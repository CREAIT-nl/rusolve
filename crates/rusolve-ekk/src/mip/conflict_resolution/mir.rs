//! Translated from SCIP 11.0.0 (Apache License 2.0, Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)). See THIRD_PARTY_NOTICES.md.
//!
//! SCIP source: conflict_resolution.c::MirReduction (lines 1375-1507);
//!   cuts.c::SCIPaggrRowAddCustomCons (lines 3154-3200), SCIPcalcMIR (lines 7934-8262),
//!   cutsTransformMIR (lines 5554-5828), doMIRBoundSubstitution (lines 5020-5100),
//!   determineBestBounds (lines 4814-5016), findMIRBestLb/findMIRBestUb (lines 4636-4810),
//!   varSection (lines 4135-4166), cutsRoundMIR (lines 6721-6938),
//!   cutsSubstituteMIR (lines 7488-7600, no-op: no LP rows), removeZerosQuad (lines 833-913),
//!   calcEfficacyDenseStorageQuad (lines 661-727); dbldblarith.h::SCIPdbldblEpsFloor2 (lines 337-354, 377).
//!
//! Single-row port: `allowlocal = FALSE`, `USEVBDS = FALSE`
//! (`vartypeusevbds = 0`), `POSTPROCESS = FALSE`, `FIXINTEGRALRHS = FALSE`,
//! `boundsfortrans = NULL`, `ignoresol = FALSE`, no LP rows in the aggregation. The branches
//! that are constant-false under these settings (local bounds, variable bounds, the
//! `fixintegralrhs` block, `postprocessCutQuad`, `boundsfortrans != NULL`) are not ported.

use crate::mip::conflict_resolution::{
    S2_BOUNDSWITCH, S2_EPSILON, S2_MAXCMIRSCALE, S2_MAXFRAC, S2_MINEFFICACY, S2_MINFRAC,
    S2_SUMEPSILON, S2ConflictRow, S2Context, compute_slack, remove_zero_vars, s2_is_binary,
    s2_is_eq, s2_is_gt, s2_is_infinity, s2_is_le, s2_is_lt, s2_is_positive, s2_lb_at, s2_ub_at,
};
use crate::prelude::*;
use crate::util::highs_cdouble::HighsCDouble;
use crate::util::highs_utils::{cpp_max, cpp_min};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// cuts.c:4097 NSECTIONS (cont, cont impl, int impl, bin impl, int, bin).
const NSECTIONS: usize = 6;

/// dbldblarith.h:42 QUAD_EPSILON (quad precision enabled).
const QUAD_EPSILON: f64 = 1e-12;

/// cuts.c:200 NONZERO(x): `COPYSIGN(1e-100, (x)) + (x)`.
fn nonzero(x: f64) -> f64 {
    1e-100_f64.copysign(x) + x
}

/// SCIP 11.0.0 (32dd425): dbldblarith.h::SCIPdbldblEpsFloor2 (line 377) with
/// SCIPdbldblFloor2 (lines 337-354). `SCIPquadprecEpsFloorQ(r, a, eps)`.
fn quad_eps_floor(a: HighsCDouble, eps: f64) -> HighsCDouble {
    // SCIPdbldblEpsFloor2: SCIPdbldblFloor2(rhi, rlo, ahi, (alo) + (eps))
    let ahi = a.hi;
    let alo = a.lo + eps;
    // line 340: __tmp_dbldbl_floor = floor((ahi) + (alo))
    let tmp_floor = (ahi + alo).floor();
    // line 341: SCIPdbldblSum21(rhi, rlo, ahi, alo, -__tmp_dbldbl_floor)
    let r = HighsCDouble { hi: ahi, lo: alo } + (-tmp_floor);
    // line 342
    if (r.hi - 1.0) + r.lo < 0.0 && r.hi + r.lo >= 0.0 {
        // lines 345-346: floor in double precision was fine
        HighsCDouble::from(tmp_floor)
    } else {
        // line 351: __tmp2_dbldbl_floor = floor((rhi) + (rlo))
        let tmp2_floor = (r.hi + r.lo).floor();
        // line 352: SCIPdbldblSum(rhi, rlo, __tmp_dbldbl_floor, __tmp2_dbldbl_floor)
        let mut res = HighsCDouble::from(0.0);
        HighsCDouble::two_sum(&mut res.hi, &mut res.lo, tmp_floor, tmp2_floor);
        res
    }
}

// ---------------------------------------------------------------------------
// MirReduction
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): conflict_resolution.c::MirReduction (lines 1375-1507)
///
/// Applies a MIR reduction to the reason row. Returns true if SCIPcalcMIR succeeded and the
/// row was replaced by the cut (the C function returns void; callers may ignore the result).
pub fn mir_reduction<C: S2Context>(
    ctx: &C,
    row: &mut S2ConflictRow,
    pos: HighsInt,
    residx: HighsInt,
    divisor: f64,
) -> bool {
    // lines 1399-1404
    debug_assert!(row.nnz() > 0);
    // line 1406
    debug_assert!(s2_is_gt(divisor, 0.0));

    // lines 1408-1409
    if !ctx.is_integral(residx) {
        return false;
    }

    let nvars = ctx.num_col();

    // lines 1420, 1428-1430: reference solution, all variables set to zero
    let mut refsol = vec![0.0f64; nvars];

    // lines 1423-1426
    let rownnz = row.nnz();
    let mut rowinds: Vec<HighsInt> = Vec::with_capacity(rownnz);
    let mut rowvals: Vec<f64> = Vec::with_capacity(rownnz);

    // lines 1433-1471: create a sparse row of the form ax <= b from the reason
    let rowrhs = -row.lhs;
    for i in 0..rownnz {
        let idx = row.inds[i];
        let coef = row.vals[idx as usize];
        rowinds.push(idx);
        rowvals.push(-coef);
        // We set the solution to the values used in propagation
        if idx == residx {
            let mut fracval: f64;
            let bnd: f64;
            // line 1447: computeSlack(set, vars, reasonrow, currbdchginfo, NULL, NULL)
            compute_slack(ctx, row, pos, None);
            if coef > 0.0 {
                // lines 1450-1452
                fracval = s2_ub_at(ctx, idx, pos, false) - row.slack / coef;
                bnd = ctx.global_ub(idx);
                fracval = cpp_min(fracval, bnd);
            } else {
                // lines 1456-1459
                debug_assert!(coef < 0.0);
                fracval = s2_lb_at(ctx, idx, pos, false) - row.slack / coef;
                bnd = ctx.global_lb(idx);
                fracval = cpp_max(fracval, bnd);
            }
            // line 1461
            refsol[idx as usize] = fracval;
        } else if coef > 0.0 {
            // line 1466
            refsol[idx as usize] = s2_ub_at(ctx, idx, pos, false);
        } else {
            // line 1468
            refsol[idx as usize] = s2_lb_at(ctx, idx, pos, false);
        }
    }

    // line 1476: SCIPaggrRowAddCustomCons(scip, aggrrow, rowinds, rowvals, rownnz, rowrhs,
    //            1.0 / divisor, 1, FALSE) into the fresh aggregation row (cuts.c:3154-3200)
    let weight = 1.0 / divisor;
    let mut aggrinds: Vec<HighsInt> = Vec::with_capacity(rownnz);
    let mut aggrvals: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); nvars];
    let mut aggrrhs = HighsCDouble::from(0.0);
    // cuts.c:3169
    debug_assert!(weight >= 0.0);
    // cuts.c:3179-3180: SCIPquadprecProdDD(quadprod, weight, rhs); SumQQ(rhs, rhs, quadprod)
    let quadprod = HighsCDouble::from(weight) * rowrhs;
    aggrrhs += quadprod;
    // cuts.c:3183-3198
    for i in 0..rownnz {
        let probindex = rowinds[i];
        let mut val = aggrvals[probindex as usize];
        if val.hi == 0.0 {
            aggrinds.push(probindex);
        }
        // cuts.c:3191-3192: SCIPquadprecProdDD(quadprod, vals[i], weight); SumQQ(val, val, quadprod)
        let quadprod = HighsCDouble::from(rowvals[i]) * weight;
        val += quadprod;
        // cuts.c:3195: QUAD_HI(val) = NONZERO(QUAD_HI(val))
        val.hi = nonzero(val.hi);
        debug_assert!(val.hi != 0.0);
        aggrvals[probindex as usize] = val;
    }

    // lines 1483-1484: SCIPcalcMIR(..., POSTPROCESS, BOUNDSWITCH, USEVBDS, FALSE, FIXINTEGRALRHS,
    //                  NULL, NULL, MINFRAC, MAXFRAC, 1.0, aggrrow, ...)
    let result = calc_mir_single_row(
        ctx,
        &aggrinds,
        &aggrvals,
        aggrrhs,
        &refsol,
        1.0,
        S2_BOUNDSWITCH,
        S2_MINFRAC,
        S2_MAXFRAC,
    );

    // lines 1486-1497
    let success = result.is_some();
    if let Some((cutinds, cutcoefs, cutrhs)) = result {
        // line 1488: assert(!cutislocal)
        // line 1489: conflictRowClear(blkmem, reasonrow, nvars)
        row.clear(nvars);
        // lines 1490-1496
        row.lhs = -cutrhs;
        for i in 0..cutinds.len() {
            row.inds.push(cutinds[i]);
            row.vals[cutinds[i] as usize] = -cutcoefs[i];
        }
    }

    // line 1500: remove variables with zero coefficient (always, as in the C)
    remove_zero_vars(row);

    success
}

// ---------------------------------------------------------------------------
// SCIPcalcMIR
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): cuts.c::SCIPcalcMIR (lines 7934-8262), for one aggregation row with
/// `postprocess = FALSE`, `vartypeusevbds = 0`, `allowlocal = FALSE`, `fixintegralrhs = FALSE`,
/// `boundsfortrans = NULL`, `cutefficacy = NULL`, no LP rows.
///
/// `aggrinds` / `aggrvals` (dense, by column) / `aggrrhs` are the aggregation row. Returns
/// `Some((cutinds, cutcoefs, cutrhs))` when `*success` is TRUE.
#[allow(clippy::too_many_arguments)]
fn calc_mir_single_row<C: S2Context>(
    ctx: &C,
    aggrinds: &[HighsInt],
    aggrvals: &[HighsCDouble],
    aggrrhs: HighsCDouble,
    refsol: &[f64],
    scale: f64,
    boundswitch: f64,
    minfrac: f64,
    maxfrac: f64,
) -> Option<(Vec<HighsInt>, Vec<f64>, f64)> {
    // line 7975
    debug_assert!(s2_is_positive(scale));

    let nvars = ctx.num_col();

    // line 7996-7997
    let nnz = aggrinds.len();
    let mut totalnnz = nnz;

    // lines 8000-8010: initialize sections (usevbds[l] = 0 since vartypeusevbds = 0)
    let mut secindices: Vec<Vec<HighsInt>> =
        (0..NSECTIONS).map(|_| Vec::with_capacity(nnz)).collect();
    let mut secnnz: [usize; NSECTIONS] = [0; NSECTIONS];
    let mut isenfint: [bool; NSECTIONS] = [false; NSECTIONS];
    let mut isimplint: [bool; NSECTIONS] = [false; NSECTIONS];
    for l in 0..NSECTIONS {
        isenfint[l] = l >= 2;
        isimplint[l] = (1..=3).contains(&l);
    }

    // lines 8022-8023: dense QUAD cutcoefs (clean buffer) and cutinds of size nvars
    let mut cutcoefs: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); nvars];
    let mut cutinds: Vec<HighsInt> = vec![0; nvars];

    // line 8025: SCIPquadprecProdQD(data->cutrhs, aggrrow->rhs, scale)
    let mut cutrhs = aggrrhs * scale;

    // lines 8027-8054
    if nnz > 0 {
        // line 8030
        cutinds[..nnz].copy_from_slice(aggrinds);

        for l in 0..nnz {
            let m = aggrinds[l] as usize;
            // line 8037: QUAD_ARRAY_LOAD(coef, aggrrow->vals, m)
            let mut coef = aggrvals[m];
            // line 8039: SCIPquadprecProdQD(coef, coef, scale)
            coef = coef * scale;
            // line 8041
            cutcoefs[m] = coef;
            debug_assert!(coef.hi != 0.0);
        }

        // line 8047: SCIPsortDownInt(data->cutinds, nnz) (indices are unique)
        cutinds[..nnz].sort_unstable_by(|a, b| b.cmp(a));
        // lines 8048-8053
        for l in 0..nnz {
            let section = var_section(ctx, cutinds[l]);
            secindices[section].push(cutinds[l]);
            secnnz[section] += 1;
        }
    }

    // line 8059-8060
    let mut ncutinds: usize = 0;
    let mut tmpislocal = false; // aggrrow->local: the custom constraint is added with local = FALSE

    // lines 8063-8064
    let mut varsign: Vec<i32> = vec![0; nvars];
    let mut boundtype: Vec<i32> = vec![0; nvars];

    // lines 8081-8094
    if totalnnz > 0 {
        let mut freevariable = false;
        let mut localbdsused = false;
        cuts_transform_mir(
            ctx,
            &mut cutcoefs,
            &mut cutrhs,
            &mut cutinds,
            &mut ncutinds,
            &mut totalnnz,
            &mut secindices,
            &mut secnnz,
            refsol,
            boundswitch,
            &mut varsign,
            &mut boundtype,
            &mut freevariable,
            &mut localbdsused,
        );
        // line 8086: assert(allowlocal || !localbdsused)
        debug_assert!(!localbdsused);
        // line 8087
        tmpislocal = tmpislocal || localbdsused;

        // lines 8089-8090
        if freevariable {
            return None;
        }
    }

    // line 8124: SCIPquadprecEpsFloorQ(downrhs, data->cutrhs, SCIPepsilon(scip))
    let downrhs = quad_eps_floor(cutrhs, S2_EPSILON);
    // line 8125: SCIPquadprecSumQQ(f0, data->cutrhs, -downrhs)
    let mut f0 = cutrhs + (-downrhs);
    // line 8126
    debug_assert!(f64::from(f0) >= -S2_EPSILON && f64::from(f0) < 1.0);

    // lines 8128-8129
    if f64::from(f0) < minfrac || f64::from(f0) > maxfrac {
        return None;
    }

    // lines 8134-8135
    if scale.abs() / (1.0 - f64::from(f0)) > S2_MAXCMIRSCALE {
        return None;
    }

    // line 8138: SCIPquadprecSumDD(f0, QUAD_HI(f0), QUAD_LO(f0))
    f0.renormalize();

    // line 8140: QUAD_ASSIGN_Q(data->cutrhs, downrhs)
    cutrhs = downrhs;

    // lines 8142-8148
    if totalnnz > 0 {
        cuts_round_mir(
            ctx,
            &mut cutcoefs,
            &mut cutrhs,
            &mut cutinds,
            &mut ncutinds,
            &mut totalnnz,
            &mut secindices,
            &mut secnnz,
            &isenfint,
            &isimplint,
            &varsign,
            &boundtype,
            f0,
        );
    }

    // lines 8164-8166: cutsSubstituteMIR with nrowinds = 0 only computes onedivoneminusf0
    // (cuts.c:7519-7520) and runs an empty loop: no effect.

    // line 8180: postprocess = FALSE
    let success = !remove_zeros_quad(
        ctx,
        S2_SUMEPSILON,
        tmpislocal,
        &mut cutcoefs,
        &mut cutrhs,
        &mut cutinds,
        &mut ncutinds,
    );

    // lines 8186-8223
    if !success {
        return None;
    }

    // line 8188
    let mirefficacy =
        calc_efficacy_dense_storage_quad(refsol, &cutcoefs, f64::from(cutrhs), &cutinds, ncutinds);

    // line 8190: SCIPisEfficacious(scip, mirefficacy) = EPSP(efficacy, sepa_minefficacy)
    // (set.c:7756-7768, root = FALSE; def.h:188 EPSP(x,eps) is (x) > (eps)); cutefficacy == NULL.
    if !(mirefficacy > S2_MINEFFICACY) {
        return None;
    }

    // lines 8192-8209
    let result_cutinds: Vec<HighsInt> = cutinds[..ncutinds].to_vec();
    let result_cutrhs = f64::from(cutrhs);
    let mut result_cutcoefs: Vec<f64> = Vec::with_capacity(ncutinds);
    for i in 0..ncutinds {
        let j = result_cutinds[i] as usize;
        result_cutcoefs.push(f64::from(cutcoefs[j]));
        cutcoefs[j] = HighsCDouble::from(0.0);
    }

    Some((result_cutinds, result_cutcoefs, result_cutrhs))
}

// ---------------------------------------------------------------------------
// cutsTransformMIR
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): cuts.c::cutsTransformMIR (lines 5554-5828), with
/// `allowlocal = FALSE`, `fixintegralrhs = FALSE` (the block at 5693-5816 is not ported),
/// `ignoresol = FALSE`, `boundsfortrans = NULL`, `usevbds[s] = 0`.
#[allow(clippy::too_many_arguments)]
fn cuts_transform_mir<C: S2Context>(
    ctx: &C,
    cutcoefs: &mut [HighsCDouble],
    cutrhs: &mut HighsCDouble,
    cutinds: &mut [HighsInt],
    ncutinds: &mut usize,
    data_totalnnz: &mut usize,
    secindices: &mut [Vec<HighsInt>],
    secnnz: &mut [usize; NSECTIONS],
    refsol: &[f64],
    boundswitch: f64,
    varsign: &mut [i32],
    boundtype: &mut [i32],
    freevariable: &mut bool,
    localbdsused: &mut bool,
) {
    // line 5592
    let totalnnz = *data_totalnnz;

    // lines 5594-5595
    *freevariable = false;
    *localbdsused = false;

    // lines 5597-5603
    let allocsize = cpp_min(NSECTIONS * totalnnz, ctx.num_col());
    let mut bestlbs: Vec<f64> = vec![0.0; allocsize];
    let mut bestubs: Vec<f64> = vec![0.0; allocsize];
    let mut bestlbtypes: Vec<i32> = vec![0; allocsize];
    let mut bestubtypes: Vec<i32> = vec![0; allocsize];
    // SCIP_BOUNDTYPE: 0 = LOWER, 1 = UPPER
    let mut selectedbounds: Vec<i32> = vec![0; allocsize];

    // line 5606: transform the cut, one variable section at a time
    for s in 0..NSECTIONS {
        let cutindsstart = *ncutinds;

        // lines 5612-5662
        let mut i: usize = 0;
        while i < secnnz[s] {
            let v = secindices[s][i];

            // lines 5619-5629: due to variable bound usage, cancellation may have occurred
            let coef = cutcoefs[v as usize];
            if f64::from(coef).abs() <= QUAD_EPSILON {
                cutcoefs[v as usize] = HighsCDouble::from(0.0);
                secnnz[s] -= 1;
                *data_totalnnz -= 1;
                secindices[s][i] = secindices[s][secnnz[s]];
                continue;
            }

            // lines 5631-5636
            let cutindex = *ncutinds;
            debug_assert!(cutindex < allocsize);
            determine_best_bounds(
                ctx,
                v,
                refsol,
                boundswitch,
                &mut bestlbs[cutindex],
                &mut bestubs[cutindex],
                &mut bestlbtypes[cutindex],
                &mut bestubtypes[cutindex],
                &mut selectedbounds[cutindex],
                freevariable,
            );

            // lines 5638-5641
            cutinds[cutindex] = v;
            *ncutinds += 1;
            i += 1;

            // lines 5644-5662: free variable, terminate
            if *freevariable {
                *ncutinds = 0;
                for j in 0..NSECTIONS {
                    for k in 0..secnnz[j] {
                        cutinds[*ncutinds] = secindices[j][k];
                        *ncutinds += 1;
                    }
                }
                return;
            }
        }

        // lines 5666-5690: perform bound substitution for added variables
        for i in cutindsstart..*ncutinds {
            let v = cutinds[i];

            let bestbnd: f64 = if selectedbounds[i] == 0 {
                debug_assert!(!s2_is_infinity(-bestlbs[i]));
                // use lower bound as transformation bound: x'_j := x_j - lb_j
                boundtype[i] = bestlbtypes[i];
                varsign[i] = 1;
                bestlbs[i]
            } else {
                debug_assert!(!s2_is_infinity(bestubs[i]));
                // use upper bound as transformation bound: x'_j := ub_j - x_j
                boundtype[i] = bestubtypes[i];
                varsign[i] = -1;
                bestubs[i]
            };
            do_mir_bound_substitution(
                cutcoefs,
                cutrhs,
                varsign[i],
                boundtype[i],
                bestbnd,
                v,
                localbdsused,
            );
        }
    }

    // lines 5693-5816: if( fixintegralrhs ) -- constant FALSE (S2_FIXINTEGRALRHS), not ported.
}

/// SCIP 11.0.0 (32dd425): cuts.c::doMIRBoundSubstitution (lines 5020-5100). Only the standard
/// bound branch (`boundtype < 0`, lines 5040-5045) is reachable: USEVBDS = FALSE.
fn do_mir_bound_substitution(
    cutcoefs: &[HighsCDouble],
    cutrhs: &mut HighsCDouble,
    varsign: i32,
    btype: i32,
    boundval: f64,
    probindex: HighsInt,
    localbdsused: &mut bool,
) {
    // line 5036
    debug_assert!(!s2_is_infinity(-(varsign as f64) * boundval));

    // line 5038
    let coef = cutcoefs[probindex as usize];

    // lines 5041-5045
    debug_assert!(btype < 0);
    // SCIPquadprecProdQD(tmp, coef, boundval); SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, -tmp)
    let tmp = coef * boundval;
    *cutrhs = *cutrhs + (-tmp);
    *localbdsused = *localbdsused || (btype == -2);
}

// ---------------------------------------------------------------------------
// determineBestBounds
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): cuts.c::determineBestBounds (lines 4814-5016), automatic bound
/// selection branch (lines 4932-4984) with `boundsfortrans = NULL`, `usevbds = 0`,
/// `allowlocal = FALSE`, `ignoresol = FALSE`.
#[allow(clippy::too_many_arguments)]
fn determine_best_bounds<C: S2Context>(
    ctx: &C,
    v: HighsInt,
    refsol: &[f64],
    boundswitch: f64,
    bestlb: &mut f64,
    bestub: &mut f64,
    bestlbtype: &mut i32,
    bestubtype: &mut i32,
    selectedbound: &mut i32,
    freevariable: &mut bool,
) {
    // lines 4939-4942
    let simplelb = find_mir_best_lb(ctx, v, bestlb, bestlbtype);
    let simpleub = find_mir_best_ub(ctx, v, bestub, bestubtype);

    // lines 4945-4952: check, if variable is free variable
    if s2_is_infinity(-*bestlb) && s2_is_infinity(*bestub) {
        *freevariable = true;
        return;
    }

    // lines 4954-4984: !ignoresol, select transformation bound
    let varsol = refsol[v as usize];

    if s2_is_infinity(*bestub) {
        *selectedbound = 0;
    } else if s2_is_infinity(-*bestlb) {
        *selectedbound = 1;
    } else if s2_is_lt(
        varsol,
        (1.0 - boundswitch) * (*bestlb) + boundswitch * (*bestub),
    ) {
        *selectedbound = 0;
    } else if s2_is_gt(
        varsol,
        (1.0 - boundswitch) * (*bestlb) + boundswitch * (*bestub),
    ) {
        *selectedbound = 1;
    } else if *bestlbtype == -1 {
        *selectedbound = 0;
    } else if *bestubtype == -1 {
        *selectedbound = 1;
    } else if (*bestlbtype >= 0 || *bestubtype >= 0)
        && !s2_is_eq(*bestlb - simplelb, simpleub - *bestub)
    {
        if *bestlb - simplelb > simpleub - *bestub {
            *selectedbound = 0;
        } else {
            *selectedbound = 1;
        }
    } else if *bestlbtype >= 0 {
        *selectedbound = 0;
    } else if *bestubtype >= 0 {
        *selectedbound = 1;
    } else {
        *selectedbound = 0;
    }
    // lines 4985-5011: ignoresol branch, constant FALSE, not ported.
}

// ---------------------------------------------------------------------------
// findMIRBestLb / findMIRBestUb
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): cuts.c::findMIRBestLb (lines 4636-4720) with `allowlocal = FALSE`
/// and `usevbds = 0` (both blocks constant-false). Returns `*simplebound`.
fn find_mir_best_lb<C: S2Context>(
    ctx: &C,
    v: HighsInt,
    bestlb: &mut f64,
    bestlbtype: &mut i32,
) -> f64 {
    *bestlb = ctx.global_lb(v);
    *bestlbtype = -1;
    // *simplebound = *bestlb
    *bestlb
}

/// SCIP 11.0.0 (32dd425): cuts.c::findMIRBestUb (lines 4724-4810) with `allowlocal = FALSE`
/// and `usevbds = 0` (both blocks constant-false). Returns `*simplebound`.
fn find_mir_best_ub<C: S2Context>(
    ctx: &C,
    v: HighsInt,
    bestub: &mut f64,
    bestubtype: &mut i32,
) -> f64 {
    *bestub = ctx.global_ub(v);
    *bestubtype = -1;
    // *simplebound = *bestub
    *bestub
}

// ---------------------------------------------------------------------------
// cutsRoundMIR
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): cuts.c::cutsRoundMIR (lines 6721-6938). Only standard global bounds
/// are used in the transformation (`type == -1`), so the local-bound and variable-bound branches
/// are unreachable.
#[allow(clippy::too_many_arguments)]
fn cuts_round_mir<C: S2Context>(
    ctx: &C,
    cutcoefs: &mut [HighsCDouble],
    cutrhs: &mut HighsCDouble,
    cutinds: &mut [HighsInt],
    ncutinds: &mut usize,
    totalnnz: &mut usize,
    secindices: &mut [Vec<HighsInt>],
    secnnz: &mut [usize; NSECTIONS],
    isenfint: &[bool; NSECTIONS],
    isimplint: &[bool; NSECTIONS],
    varsign: &[i32],
    boundtype: &[i32],
    f0: HighsCDouble,
) {
    // line 6738
    debug_assert!(0.0 < f64::from(f0) && f64::from(f0) < 1.0);

    // line 6740: SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0)
    let mut onedivoneminusf0 = (-f0) + 1.0;
    // line 6741: SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0)
    onedivoneminusf0 = 1.0 / onedivoneminusf0;

    // line 6746
    let mut cutindex: isize = *ncutinds as isize - 1;
    // line 6747
    for s in (0..NSECTIONS).rev() {
        let nnz = secnnz[s];
        let enfintegral = isenfint[s];
        let implintegral = isimplint[s];

        // line 6756
        for i in (0..nnz).rev() {
            let v = secindices[s][i];
            debug_assert!(0 <= v && (v as usize) < ctx.num_col());
            debug_assert!(cutinds[cutindex as usize] == v);

            let sign = varsign[cutindex as usize];
            debug_assert!(sign == 1 || sign == -1);
            let vtype = boundtype[cutindex as usize];

            cutindex -= 1;

            // line 6779
            let mut aj = cutcoefs[v as usize];
            let mut cutaj: HighsCDouble;

            if enfintegral || implintegral {
                // line 6787: QUAD_SCALE(aj, sign)
                aj.hi *= sign as f64;
                aj.lo *= sign as f64;

                // line 6789: SCIPquadprecEpsFloorQ(downaj, aj, SCIPepsilon(scip))
                let downaj = quad_eps_floor(aj, S2_EPSILON);
                // line 6790: SCIPquadprecSumQQ(fj, aj, -downaj)
                let fj = aj + (-downaj);
                debug_assert!(f64::from(fj) >= -S2_EPSILON && f64::from(fj) < 1.0);

                // line 6793
                if s2_is_le(f64::from(fj), f64::from(f0)) {
                    cutaj = downaj;
                } else {
                    // lines 6799-6801
                    let mut tmp = fj + (-f0);
                    tmp = tmp * onedivoneminusf0;
                    cutaj = tmp + downaj;
                }
                // line 6803: QUAD_SCALE(cutaj, sign)
                cutaj.hi *= sign as f64;
                cutaj.lo *= sign as f64;
            } else {
                // lines 6808-6813: variable is continuous
                if f64::from(aj) * (sign as f64) >= 0.0 {
                    cutaj = HighsCDouble::from(0.0);
                } else {
                    cutaj = onedivoneminusf0 * aj;
                }
            }

            // lines 6817-6825: remove coefficient from cut if it becomes zero
            if f64::from(cutaj).abs() <= QUAD_EPSILON {
                cutcoefs[v as usize] = HighsCDouble::from(0.0);
                *totalnnz -= 1;
                secnnz[s] -= 1;
                secindices[s][i] = secindices[s][secnnz[s]];
                continue;
            }

            // line 6828
            cutcoefs[v as usize] = cutaj;

            // lines 6831-6866: undo bound transformations (standard global bound, type == -1)
            debug_assert!(vtype == -1);
            if sign == 1 {
                // lines 6841-6842
                let tmp = cutaj * ctx.global_lb(v);
                *cutrhs = *cutrhs + tmp;
            } else {
                // lines 6857-6858
                let tmp = cutaj * ctx.global_ub(v);
                *cutrhs = *cutrhs + tmp;
            }
        }
    }

    // lines 6926-6936
    *ncutinds = 0;
    for s in 0..NSECTIONS {
        for i in 0..secnnz[s] {
            cutinds[*ncutinds] = secindices[s][i];
            *ncutinds += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// removeZerosQuad
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): cuts.c::removeZerosQuad (lines 833-913). Returns TRUE if the cut
/// became redundant.
fn remove_zeros_quad<C: S2Context>(
    ctx: &C,
    minval: f64,
    cutislocal: bool,
    cutcoefs: &mut [HighsCDouble],
    cutrhs: &mut HighsCDouble,
    cutinds: &mut [HighsInt],
    cutnnz: &mut usize,
) -> bool {
    // MirReduction only produces global cuts (allowlocal = FALSE), so the local-bound branch
    // (lines 859-863) is unreachable.
    debug_assert!(!cutislocal);

    let mut i: usize = 0;
    while i < *cutnnz {
        let v = cutinds[i];
        let mut val = cutcoefs[v as usize];

        // lines 866-867
        let lb = ctx.global_lb(v);
        let ub = ctx.global_ub(v);

        // lines 870-873
        let isfixed = !(s2_is_infinity(-lb) || s2_is_infinity(ub)) && s2_is_eq(ub, lb);

        // line 875: isfixed || EPSZ(QUAD_TO_DBL(val), minval)
        if isfixed || f64::from(val).abs() <= minval {
            if f64::from(val).abs() > QUAD_EPSILON {
                // adjust right hand side with max contribution
                if f64::from(val) < 0.0 {
                    if s2_is_infinity(ub) {
                        return true;
                    }
                    val = val * ub;
                    *cutrhs = *cutrhs + (-val);
                } else {
                    if s2_is_infinity(-lb) {
                        return true;
                    }
                    val = val * lb;
                    *cutrhs = *cutrhs + (-val);
                }
            }

            // lines 902-907
            cutcoefs[v as usize] = HighsCDouble::from(0.0);
            *cutnnz -= 1;
            cutinds[i] = cutinds[*cutnnz];
        } else {
            i += 1;
        }
    }

    // lines 913-914: relax rhs to 0, if it's very close to 0
    if f64::from(*cutrhs) < 0.0 && f64::from(*cutrhs) >= -S2_EPSILON {
        *cutrhs = HighsCDouble::from(0.0);
    }

    false
}

// ---------------------------------------------------------------------------
// calcEfficacyDenseStorageQuad
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): cuts.c::calcEfficacyDenseStorageQuad (lines 661-727), norm 'e'
/// (set.c:463).
fn calc_efficacy_dense_storage_quad(
    refsol: &[f64],
    cutcoefs: &[HighsCDouble],
    cutrhs: f64,
    cutinds: &[HighsInt],
    cutnnz: usize,
) -> f64 {
    let mut norm: f64 = 0.0;
    let mut activity: f64 = 0.0;

    // lines 686-692
    for i in 0..cutnnz {
        let coef = cutcoefs[cutinds[i] as usize];
        activity += f64::from(coef) * refsol[cutinds[i] as usize];
        norm += f64::from(coef) * f64::from(coef);
    }
    norm = norm.sqrt();

    // line 727: (activity - cutrhs) / MAX(1e-6, norm)
    (activity - cutrhs) / cpp_max(1e-6, norm)
}

// ---------------------------------------------------------------------------
// varSection
// ---------------------------------------------------------------------------

/// SCIP 11.0.0 (32dd425): cuts.c::varSection (lines 4135-4166).
///
/// SCIP derives the section from the problem index, because its variables are sorted by type
/// (binary, integer, implied integer, continuous). HiGHS columns are not sorted, and the
/// S2Context has no implied-integer information, so the section is taken from the column type:
/// 0 = continuous, 4 = integer, 5 = binary (SCIPvarIsBinary, mod.rs s2_is_binary). The
/// rounding treats sections 2-5 alike (enforced integral), so implied integers are rounded as
/// in SCIP; only the section order of the cut indices can differ.
fn var_section<C: S2Context>(ctx: &C, v: HighsInt) -> usize {
    if !ctx.is_integral(v) {
        return 0;
    }
    if s2_is_binary(ctx, v) {
        return 5;
    }
    4
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mip::conflict_resolution::{S2ConflictRow, S2TestChg, S2TestCtx};

    /// s2_mir_basic: 2x0 + 3x1 >= 5, x0 integer [0,3], x1 binary.
    /// Change x1 ub 1->0 at pos 0, x0 lb 0->3 at pos 1.
    /// Reason row resolved on x0 (pos 1).
    /// Expected: x0 + 2x1 >= 3, lhs 3.0, vals[0]=1, vals[1]=2.
    #[test]
    fn s2_mir_basic() {
        // x0: integer, global [0, 3]
        // x1: binary, global [0, 1]
        // Row: 2*x0 + 3*x1 >= 5
        // Changes: pos 0: x1 ub 1->0; pos 1: x0 lb 0->3
        // Resolve at pos 1 (x0 lb change)

        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0],
            gub: vec![3.0, 1.0],
            integral: vec![true, true],
            num_integral: 2,
            stack: vec![
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 3.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
            ],
        };

        let mut row = S2ConflictRow {
            vals: vec![2.0, 3.0],
            inds: vec![0, 1],
            lhs: 5.0,
            slack: 0.0,
            coefquotient: 0.0,
        };

        let residx: HighsInt = 0;
        let pos: HighsInt = 1;
        let divisor: f64 = 2.0; // coefinrow of x0 (conflict_resolution.c:2292 passes coefinrow)

        let changed = mir_reduction(&ctx, &mut row, pos, residx, divisor);

        // The row should be modified
        // Expected: x0 + 2x1 >= 3
        // lhs = 3.0, vals[0] = 1.0, vals[1] = 2.0
        assert!(changed, "mir_reduction should have modified the row");
        assert_eq!(row.lhs, 3.0, "lhs should be 3.0");
        assert_eq!(row.vals[0], 1.0, "vals[0] should be 1.0");
        assert_eq!(row.vals[1], 2.0, "vals[1] should be 2.0");
    }

    /// s2_mir_f0_zero: 2x0 + 3x1 >= 4, x0 integer [0,3], x1 binary.
    /// Change x1 ub 1->0 at pos 0, x0 lb 0->2 at pos 1.
    /// refsol x0 = 2, scaled gives x0' - 1.5x1' <= 1, f0 = 0.
    /// Row unchanged.
    #[test]
    fn s2_mir_f0_zero() {
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0],
            gub: vec![3.0, 1.0],
            integral: vec![true, true],
            num_integral: 2,
            stack: vec![
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 2.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
            ],
        };

        let mut row = S2ConflictRow {
            vals: vec![2.0, 3.0],
            inds: vec![0, 1],
            lhs: 4.0,
            slack: 0.0,
            coefquotient: 0.0,
        };

        let residx: HighsInt = 0;
        let pos: HighsInt = 1;
        let divisor: f64 = 2.0; // coefinrow of x0 (conflict_resolution.c:2292 passes coefinrow)

        let changed = mir_reduction(&ctx, &mut row, pos, residx, divisor);

        // f0 = 0, which is < MINFRAC (0.05), so the row is unchanged
        assert!(
            !changed,
            "mir_reduction should not have modified the row (f0=0)"
        );
        assert_eq!(row.lhs, 4.0, "lhs should remain 4.0");
        assert_eq!(row.vals[0], 2.0, "vals[0] should remain 2.0");
        assert_eq!(row.vals[1], 3.0, "vals[1] should remain 3.0");
    }

    /// s2_mir_continuous_resvar: the resolved column is continuous,
    /// so the row is unchanged.
    #[test]
    fn s2_mir_continuous_resvar() {
        let ctx = S2TestCtx {
            glb: vec![0.0, 0.0],
            gub: vec![3.0, 1.0],
            integral: vec![false, true],
            num_integral: 1,
            stack: vec![
                S2TestChg {
                    col: 1,
                    is_upper: true,
                    newbound: 0.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
                S2TestChg {
                    col: 0,
                    is_upper: false,
                    newbound: 3.0,
                    depth: 0,
                    branching: false,
                    reason: None,
                },
            ],
        };

        let mut row = S2ConflictRow {
            vals: vec![2.0, 3.0],
            inds: vec![0, 1],
            lhs: 5.0,
            slack: 0.0,
            coefquotient: 0.0,
        };

        let residx: HighsInt = 0; // continuous
        let pos: HighsInt = 1;
        let divisor: f64 = 2.0; // coefinrow of x0 (conflict_resolution.c:2292 passes coefinrow)

        let changed = mir_reduction(&ctx, &mut row, pos, residx, divisor);

        // residx is not integral, so the row is unchanged
        assert!(
            !changed,
            "mir_reduction should not modify row for continuous resvar"
        );
        assert_eq!(row.lhs, 5.0, "lhs should remain 5.0");
        assert_eq!(row.vals[0], 2.0, "vals[0] should remain 2.0");
        assert_eq!(row.vals[1], 3.0, "vals[1] should remain 3.0");
    }
}
