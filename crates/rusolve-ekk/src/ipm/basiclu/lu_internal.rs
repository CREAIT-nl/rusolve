//! HiGHS source: ipm/basiclu/lu_internal.h (struct lu), lu_internal.c
//!
//! Copyright (C) 2016-2018  ERGO-Code (BASICLU source)
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::*;

/// rusolve BASICLU time limit (not in the HiGHS reference): lu_internal.c::luTime (= basiclu_wallclock): wall clock in seconds.
pub fn lu_time() -> f64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs_f64())
        .unwrap_or(0.0)
}

// lu_internal.c: private entries in xstore
pub const BASICLU_TASK: usize = 256;
pub const BASICLU_FTCOLUMN_IN: usize = 257;
pub const BASICLU_FTCOLUMN_OUT: usize = 258;
pub const BASICLU_PIVOT_ROW: usize = 259;
pub const BASICLU_PIVOT_COL: usize = 260;
pub const BASICLU_RANKDEF: usize = 261;
pub const BASICLU_MIN_COLNZ: usize = 262;
pub const BASICLU_MIN_ROWNZ: usize = 263;
pub const BASICLU_MARKER: usize = 266;
pub const BASICLU_UPDATE_COST_NUMER: usize = 267;
pub const BASICLU_UPDATE_COST_DENOM: usize = 268;
pub const BASICLU_PIVOTLEN: usize = 269;

/// HiGHS: basiclu/lu_internal.h::struct lu (lines 20-139)
///
/// Scalars are held by value (copied from `xstore` by [`lu_load`], written
/// back by [`lu_save`]). Array members are slices into the caller's
/// storage. `istore[1..]` is carved in the C partition order of `lu_load`;
/// `xstore[0..512]` is `xhead` (the scalar block `lu_save` writes) and
/// `xstore[512..]` holds `work0`, `work1`, `col_pivot`, `row_pivot`.
///
/// C second-list aliases (same memory, used in solve/update instead of
/// factorize) are not fields. Translate them as:
///
/// | C alias      | Rust                         |
/// |--------------|------------------------------|
/// | `pivotcol`   | `colcount_flink`             |
/// | `pivotrow`   | `colcount_blink`             |
/// | `Rbegin`     | `rowcount_flink`             |
/// | `eta_row`    | `rowcount_flink[m+1..]`      |
/// | `iwork1`     | `rowcount_blink`             |
/// | `Lbegin`     | `wbegin[m+1..]`              |
/// | `Ltbegin`    | `wend[m+1..]`                |
/// | `Ltbegin_p`  | `wflink[m+1..]`              |
/// | `p`          | `wblink[m+1..]`              |
/// | `pmap`       | `pinv`                       |
/// | `qmap`       | `qinv`                       |
/// | `marked`     | `iwork0`                     |
///
/// When one function needs a region and its alias at once (`Wbegin` and
/// `Lbegin`), split the region at the C offset with `split_at_mut`.
pub struct Lu<'a> {
    // user parameters, not modified
    pub lmem: LuInt,
    pub umem: LuInt,
    pub wmem: LuInt,
    pub droptol: f64,
    pub abstol: f64,
    pub reltol: f64,
    pub nzbias: LuInt,
    pub maxsearch: LuInt,
    pub pad: LuInt,
    pub stretch: f64,
    pub compress_thres: f64,
    pub sparse_thres: f64,
    pub search_rows: LuInt,
    /// rusolve BASICLU time limit (not in the HiGHS reference): lu_internal.h::basiclu_time_start / basiclu_time_limit
    pub basiclu_time_start: f64,
    pub basiclu_time_limit: f64,

    // user readable
    pub m: LuInt,
    pub addmem_l: LuInt,
    pub addmem_u: LuInt,
    pub addmem_w: LuInt,

    pub nupdate: LuInt,
    pub nforrest: LuInt,
    pub nfactorize: LuInt,
    pub nupdate_total: LuInt,
    pub nforrest_total: LuInt,
    pub nsymperm_total: LuInt,
    pub lnz: LuInt,
    pub unz: LuInt,
    pub rnz: LuInt,
    pub min_pivot: f64,
    pub max_pivot: f64,
    pub max_eta: f64,
    pub update_cost_numer: f64,
    pub update_cost_denom: f64,
    pub time_factorize: f64,
    pub time_solve: f64,
    pub time_update: f64,
    pub time_factorize_total: f64,
    pub time_solve_total: f64,
    pub time_update_total: f64,
    pub lflops: LuInt,
    pub uflops: LuInt,
    pub rflops: LuInt,
    pub condest_l: f64,
    pub condest_u: f64,
    pub norm_l: f64,
    pub norm_u: f64,
    pub normest_linv: f64,
    pub normest_uinv: f64,
    pub onenorm: f64,
    pub infnorm: f64,
    pub residual_test: f64,

    pub matrix_nz: LuInt,
    pub rank: LuInt,
    pub bump_size: LuInt,
    pub bump_nz: LuInt,
    pub nsearch_pivot: LuInt,
    pub nexpand: LuInt,
    pub ngarbage: LuInt,
    pub factor_flops: LuInt,
    pub time_singletons: f64,
    pub time_search_pivot: f64,
    pub time_elim_pivot: f64,

    pub pivot_error: f64,

    // private
    pub task: LuInt,
    pub pivot_row: LuInt,
    pub pivot_col: LuInt,
    pub ftran_for_update: LuInt,
    pub btran_for_update: LuInt,
    pub marker: LuInt,
    pub pivotlen: LuInt,
    pub rankdef: LuInt,
    pub min_colnz: LuInt,
    pub min_rownz: LuInt,

    // aliases to user arrays
    pub lindex: &'a mut [LuInt],
    pub uindex: &'a mut [LuInt],
    pub windex: &'a mut [LuInt],
    pub lvalue: &'a mut [f64],
    pub uvalue: &'a mut [f64],
    pub wvalue: &'a mut [f64],

    // regions of istore[1..] (C partition order)
    /// size 2m+2
    pub colcount_flink: &'a mut [LuInt],
    /// size 2m+2
    pub colcount_blink: &'a mut [LuInt],
    /// size 2m+2
    pub rowcount_flink: &'a mut [LuInt],
    /// size 2m+2
    pub rowcount_blink: &'a mut [LuInt],
    /// size 2m+1
    pub wbegin: &'a mut [LuInt],
    /// size 2m+1
    pub wend: &'a mut [LuInt],
    /// size 2m+1
    pub wflink: &'a mut [LuInt],
    /// size 2m+1
    pub wblink: &'a mut [LuInt],
    /// size m
    pub pinv: &'a mut [LuInt],
    /// size m
    pub qinv: &'a mut [LuInt],
    /// size m+1
    pub lbegin_p: &'a mut [LuInt],
    /// size m+1
    pub ubegin: &'a mut [LuInt],
    /// size m, zeroed
    pub iwork0: &'a mut [LuInt],

    // regions of xstore
    /// `xstore[0..512]`: the scalar block written by [`lu_save`]
    pub xhead: &'a mut [f64],
    /// size m workspace, zeroed
    pub work0: &'a mut [f64],
    /// size m workspace, uninitialized (also the `pstack` pun)
    pub work1: &'a mut [f64],
    /// pivot elements by column index
    pub col_pivot: &'a mut [f64],
    /// pivot elements by row index (also the `cancelled` pun in `lu_pivot_small`)
    pub row_pivot: &'a mut [f64],
}

/// HiGHS: basiclu/lu_internal.c::lu_load (lines 36-191)
///
/// PORT NOTE: C fills a caller-declared `struct lu` and returns a status;
/// here the loaded [`Lu`] is returned, or the error status
/// (`BASICLU_ERROR_invalid_store`) with nothing loaded. `xstore` values are
/// converted to integers with `as LuInt` (C implicit `double` -> `lu_int`).
pub fn lu_load<'a>(
    istore: &'a mut [LuInt],
    xstore: &'a mut [f64],
    li: &'a mut [LuInt],
    lx: &'a mut [f64],
    ui: &'a mut [LuInt],
    ux: &'a mut [f64],
    wi: &'a mut [LuInt],
    wx: &'a mut [f64],
) -> Result<Lu<'a>, LuInt> {
    // C: if (!istore || istore[0] != BASICLU_HASH || !xstore || xstore[0] != BASICLU_HASH)
    // PORT NOTE: an empty `istore` plays the C null pointer. A store shorter than the
    // 512-entry scalar block is rejected too (C would read past its end).
    if istore.is_empty()
        || istore[0] != BASICLU_HASH
        || xstore.len() < 512
        || xstore[0] != BASICLU_HASH as f64
    {
        return Err(BASICLU_ERROR_invalid_store);
    }

    // xstore[0..512] is the scalar block (xhead); xstore[512..] holds work0..row_pivot.
    let (xhead, xwork) = xstore.split_at_mut(512);
    // C implicit double -> lu_int conversion is `as LuInt`.
    let xi = |k: usize| xhead[k] as LuInt;
    let xd = |k: usize| xhead[k];

    let m: LuInt = xi(BASICLU_DIM);
    let m_us = m as usize;

    // partition istore for factorize (lu_internal.c 136-149); istore[0] is the hash
    let (_, rest) = istore.split_at_mut(1);
    let (colcount_flink, rest) = rest.split_at_mut(2 * m_us + 2);
    let (colcount_blink, rest) = rest.split_at_mut(2 * m_us + 2);
    let (rowcount_flink, rest) = rest.split_at_mut(2 * m_us + 2);
    let (rowcount_blink, rest) = rest.split_at_mut(2 * m_us + 2);
    let (wbegin, rest) = rest.split_at_mut(2 * m_us + 1);
    let (wend, rest) = rest.split_at_mut(2 * m_us + 1);
    let (wflink, rest) = rest.split_at_mut(2 * m_us + 1);
    let (wblink, rest) = rest.split_at_mut(2 * m_us + 1);
    let (pinv, rest) = rest.split_at_mut(m_us);
    let (qinv, rest) = rest.split_at_mut(m_us);
    let (lbegin_p, rest) = rest.split_at_mut(m_us + 1);
    let (ubegin, rest) = rest.split_at_mut(m_us + 1);
    let (iwork0, _) = rest.split_at_mut(m_us);

    // partition xstore for factorize and update (lu_internal.c 165-170)
    let (work0, xwork) = xwork.split_at_mut(m_us);
    let (work1, xwork) = xwork.split_at_mut(m_us);
    let (col_pivot, xwork) = xwork.split_at_mut(m_us);
    let (row_pivot, _) = xwork.split_at_mut(m_us);

    // Reset @marked if increasing @marker by four causes overflow.
    let mut marker = xi(BASICLU_MARKER);
    if marker > LU_INT_MAX - 4 {
        iwork0.fill(0); // memset(marked, 0, m*sizeof(lu_int))
        marker = 0;
    }

    // One past the final position in @Wend must hold the file size.
    // The file has 2*m lines while factorizing and m lines otherwise.
    let nupdate = xi(BASICLU_NUPDATE);
    let wmem = xi(BASICLU_MEMORYW);
    if nupdate >= 0 {
        wend[m_us] = wmem;
    } else {
        wend[2 * m_us] = wmem;
    }

    Ok(Lu {
        lmem: xi(BASICLU_MEMORYL),
        umem: xi(BASICLU_MEMORYU),
        wmem,
        droptol: xd(BASICLU_DROP_TOLERANCE),
        abstol: xd(BASICLU_ABS_PIVOT_TOLERANCE),
        reltol: xd(BASICLU_REL_PIVOT_TOLERANCE).min(1.0),
        nzbias: xi(BASICLU_BIAS_NONZEROS),
        maxsearch: xi(BASICLU_MAXN_SEARCH_PIVOT),
        pad: xi(BASICLU_PAD),
        stretch: xd(BASICLU_STRETCH),
        compress_thres: xd(BASICLU_COMPRESSION_THRESHOLD),
        sparse_thres: xd(BASICLU_SPARSE_THRESHOLD),
        search_rows: (xd(BASICLU_SEARCH_ROWS) != 0.0) as LuInt,
        basiclu_time_start: xd(BASICLU_TIME_START),
        basiclu_time_limit: xd(BASICLU_TIME_LIMIT),

        m,
        addmem_l: 0,
        addmem_u: 0,
        addmem_w: 0,

        nupdate,
        nforrest: xi(BASICLU_NFORREST),
        nfactorize: xi(BASICLU_NFACTORIZE),
        nupdate_total: xi(BASICLU_NUPDATE_TOTAL),
        nforrest_total: xi(BASICLU_NFORREST_TOTAL),
        nsymperm_total: xi(BASICLU_NSYMPERM_TOTAL),
        lnz: xi(BASICLU_LNZ),
        unz: xi(BASICLU_UNZ),
        rnz: xi(BASICLU_RNZ),
        min_pivot: xd(BASICLU_MIN_PIVOT),
        max_pivot: xd(BASICLU_MAX_PIVOT),
        max_eta: xd(BASICLU_MAX_ETA),
        update_cost_numer: xd(BASICLU_UPDATE_COST_NUMER),
        update_cost_denom: xd(BASICLU_UPDATE_COST_DENOM),
        time_factorize: xd(BASICLU_TIME_FACTORIZE),
        time_solve: xd(BASICLU_TIME_SOLVE),
        time_update: xd(BASICLU_TIME_UPDATE),
        time_factorize_total: xd(BASICLU_TIME_FACTORIZE_TOTAL),
        time_solve_total: xd(BASICLU_TIME_SOLVE_TOTAL),
        time_update_total: xd(BASICLU_TIME_UPDATE_TOTAL),
        lflops: xi(BASICLU_LFLOPS),
        uflops: xi(BASICLU_UFLOPS),
        rflops: xi(BASICLU_RFLOPS),
        condest_l: xd(BASICLU_CONDEST_L),
        condest_u: xd(BASICLU_CONDEST_U),
        norm_l: xd(BASICLU_NORM_L),
        norm_u: xd(BASICLU_NORM_U),
        normest_linv: xd(BASICLU_NORMEST_LINV),
        normest_uinv: xd(BASICLU_NORMEST_UINV),
        onenorm: xd(BASICLU_MATRIX_ONENORM),
        infnorm: xd(BASICLU_MATRIX_INFNORM),
        residual_test: xd(BASICLU_RESIDUAL_TEST),

        matrix_nz: xi(BASICLU_MATRIX_NZ),
        rank: xi(BASICLU_RANK),
        bump_size: xi(BASICLU_BUMP_SIZE),
        bump_nz: xi(BASICLU_BUMP_NZ),
        nsearch_pivot: xi(BASICLU_NSEARCH_PIVOT),
        nexpand: xi(BASICLU_NEXPAND),
        ngarbage: xi(BASICLU_NGARBAGE),
        factor_flops: xi(BASICLU_FACTOR_FLOPS),
        time_singletons: xd(BASICLU_TIME_SINGLETONS),
        time_search_pivot: xd(BASICLU_TIME_SEARCH_PIVOT),
        time_elim_pivot: xd(BASICLU_TIME_ELIM_PIVOT),

        pivot_error: xd(BASICLU_PIVOT_ERROR),

        task: xi(BASICLU_TASK),
        pivot_row: xi(BASICLU_PIVOT_ROW),
        pivot_col: xi(BASICLU_PIVOT_COL),
        ftran_for_update: xi(BASICLU_FTCOLUMN_IN),
        btran_for_update: xi(BASICLU_FTCOLUMN_OUT),
        marker,
        pivotlen: xi(BASICLU_PIVOTLEN),
        rankdef: xi(BASICLU_RANKDEF),
        min_colnz: xi(BASICLU_MIN_COLNZ),
        min_rownz: xi(BASICLU_MIN_ROWNZ),

        // aliases to user arrays (C: Lindex = Li, ...)
        lindex: li,
        uindex: ui,
        windex: wi,
        lvalue: lx,
        uvalue: ux,
        wvalue: wx,

        colcount_flink,
        colcount_blink,
        rowcount_flink,
        rowcount_blink,
        wbegin,
        wend,
        wflink,
        wblink,
        pinv,
        qinv,
        lbegin_p,
        ubegin,
        iwork0,

        xhead,
        work0,
        work1,
        col_pivot,
        row_pivot,
    })
}

/// HiGHS: basiclu/lu_internal.c::lu_save (lines 203-274)
///
/// PORT NOTE: C takes `istore, xstore`; the scalar block is `this.xhead`.
pub fn lu_save(this: &mut Lu<'_>, status: LuInt) -> LuInt {
    // user readable
    this.xhead[BASICLU_STATUS] = status as f64;
    this.xhead[BASICLU_ADD_MEMORYL] = this.addmem_l as f64;
    this.xhead[BASICLU_ADD_MEMORYU] = this.addmem_u as f64;
    this.xhead[BASICLU_ADD_MEMORYW] = this.addmem_w as f64;

    this.xhead[BASICLU_NUPDATE] = this.nupdate as f64;
    this.xhead[BASICLU_NFORREST] = this.nforrest as f64;
    this.xhead[BASICLU_NFACTORIZE] = this.nfactorize as f64;
    this.xhead[BASICLU_NUPDATE_TOTAL] = this.nupdate_total as f64;
    this.xhead[BASICLU_NFORREST_TOTAL] = this.nforrest_total as f64;
    this.xhead[BASICLU_NSYMPERM_TOTAL] = this.nsymperm_total as f64;
    this.xhead[BASICLU_LNZ] = this.lnz as f64;
    this.xhead[BASICLU_UNZ] = this.unz as f64;
    this.xhead[BASICLU_RNZ] = this.rnz as f64;
    this.xhead[BASICLU_MIN_PIVOT] = this.min_pivot;
    this.xhead[BASICLU_MAX_PIVOT] = this.max_pivot;
    this.xhead[BASICLU_MAX_ETA] = this.max_eta;
    this.xhead[BASICLU_UPDATE_COST_NUMER] = this.update_cost_numer;
    this.xhead[BASICLU_UPDATE_COST_DENOM] = this.update_cost_denom;
    this.xhead[BASICLU_UPDATE_COST] = this.update_cost_numer / this.update_cost_denom;
    this.xhead[BASICLU_TIME_FACTORIZE] = this.time_factorize;
    this.xhead[BASICLU_TIME_SOLVE] = this.time_solve;
    this.xhead[BASICLU_TIME_UPDATE] = this.time_update;
    this.xhead[BASICLU_TIME_FACTORIZE_TOTAL] = this.time_factorize_total;
    this.xhead[BASICLU_TIME_SOLVE_TOTAL] = this.time_solve_total;
    this.xhead[BASICLU_TIME_UPDATE_TOTAL] = this.time_update_total;
    this.xhead[BASICLU_LFLOPS] = this.lflops as f64;
    this.xhead[BASICLU_UFLOPS] = this.uflops as f64;
    this.xhead[BASICLU_RFLOPS] = this.rflops as f64;
    this.xhead[BASICLU_CONDEST_L] = this.condest_l;
    this.xhead[BASICLU_CONDEST_U] = this.condest_u;
    this.xhead[BASICLU_NORM_L] = this.norm_l;
    this.xhead[BASICLU_NORM_U] = this.norm_u;
    this.xhead[BASICLU_NORMEST_LINV] = this.normest_linv;
    this.xhead[BASICLU_NORMEST_UINV] = this.normest_uinv;
    this.xhead[BASICLU_MATRIX_ONENORM] = this.onenorm;
    this.xhead[BASICLU_MATRIX_INFNORM] = this.infnorm;
    this.xhead[BASICLU_RESIDUAL_TEST] = this.residual_test;

    this.xhead[BASICLU_MATRIX_NZ] = this.matrix_nz as f64;
    this.xhead[BASICLU_RANK] = this.rank as f64;
    this.xhead[BASICLU_BUMP_SIZE] = this.bump_size as f64;
    this.xhead[BASICLU_BUMP_NZ] = this.bump_nz as f64;
    this.xhead[BASICLU_NSEARCH_PIVOT] = this.nsearch_pivot as f64;
    this.xhead[BASICLU_NEXPAND] = this.nexpand as f64;
    this.xhead[BASICLU_NGARBAGE] = this.ngarbage as f64;
    this.xhead[BASICLU_FACTOR_FLOPS] = this.factor_flops as f64;
    this.xhead[BASICLU_TIME_SINGLETONS] = this.time_singletons;
    this.xhead[BASICLU_TIME_SEARCH_PIVOT] = this.time_search_pivot;
    this.xhead[BASICLU_TIME_ELIM_PIVOT] = this.time_elim_pivot;

    this.xhead[BASICLU_PIVOT_ERROR] = this.pivot_error;

    // private
    this.xhead[BASICLU_TASK] = this.task as f64;
    this.xhead[BASICLU_PIVOT_ROW] = this.pivot_row as f64;
    this.xhead[BASICLU_PIVOT_COL] = this.pivot_col as f64;
    this.xhead[BASICLU_FTCOLUMN_IN] = this.ftran_for_update as f64;
    this.xhead[BASICLU_FTCOLUMN_OUT] = this.btran_for_update as f64;
    this.xhead[BASICLU_MARKER] = this.marker as f64;
    this.xhead[BASICLU_PIVOTLEN] = this.pivotlen as f64;
    this.xhead[BASICLU_RANKDEF] = this.rankdef as f64;
    this.xhead[BASICLU_MIN_COLNZ] = this.min_colnz as f64;
    this.xhead[BASICLU_MIN_ROWNZ] = this.min_rownz as f64;

    status
}

/// HiGHS: basiclu/lu_internal.c::lu_reset (lines 283-352)
pub fn lu_reset(this: &mut Lu<'_>) {
    // user readable
    this.nupdate = -1; // invalidate factorization
    this.nforrest = 0;
    this.lnz = 0;
    this.unz = 0;
    this.rnz = 0;
    this.min_pivot = 0.0;
    this.max_pivot = 0.0;
    this.max_eta = 0.0;
    this.update_cost_numer = 0.0;
    this.update_cost_denom = 1.0;
    this.time_factorize = 0.0;
    this.time_solve = 0.0;
    this.time_update = 0.0;
    this.lflops = 0;
    this.uflops = 0;
    this.rflops = 0;
    this.condest_l = 0.0;
    this.condest_u = 0.0;
    this.norm_l = 0.0;
    this.norm_u = 0.0;
    this.normest_linv = 0.0;
    this.normest_uinv = 0.0;
    this.onenorm = 0.0;
    this.infnorm = 0.0;
    this.residual_test = 0.0;

    this.matrix_nz = 0;
    this.rank = 0;
    this.bump_size = 0;
    this.bump_nz = 0;
    this.nsearch_pivot = 0;
    this.nexpand = 0;
    this.ngarbage = 0;
    this.factor_flops = 0;
    this.time_singletons = 0.0;
    this.time_search_pivot = 0.0;
    this.time_elim_pivot = 0.0;

    this.pivot_error = 0.0;

    // private
    this.task = NO_TASK;
    this.pivot_row = -1;
    this.pivot_col = -1;
    this.ftran_for_update = -1;
    this.btran_for_update = -1;
    this.marker = 0;
    this.pivotlen = 0;
    this.rankdef = 0;
    this.min_colnz = 1;
    this.min_rownz = 1;

    let m = this.m as usize;

    // One past the final position in @Wend must hold the file size.
    // The file has 2*m lines during factorization.
    this.wend[2 * m] = this.wmem;

    // The integer workspace iwork0 must be zeroed for a new factorization.
    // The double workspace work0 actually needs only be zeroed once in the
    // initialization of xstore. However, it is easier and more consistent
    // to do that here as well.
    this.iwork0[..m].fill(0);

    // BASICLU quirk: memset(work0, 0, m * sizeof(lu_int)) clears m*4 BYTES,
    // not m doubles. So only the first m/2 doubles are zeroed, plus the low
    // (little-endian) bytes of the next one when m is odd. The bytes are
    // cleared as bit patterns, exactly as memset does.
    let nbytes = 4 * m;
    let full = nbytes / 8;
    let rem = nbytes % 8;
    this.work0[..full].fill(0.0);
    if rem > 0 {
        let bits = this.work0[full].to_bits();
        this.work0[full] = f64::from_bits(bits & !((1u64 << (8 * rem)) - 1));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A valid BASICLU store with dimension `m`; factor arrays of length `nf`.
    fn fresh_store(m: LuInt) -> (Vec<LuInt>, Vec<f64>) {
        let mut istore = vec![0; 1024];
        let mut xstore = vec![0.0; 1024];
        istore[0] = BASICLU_HASH;
        xstore[0] = BASICLU_HASH as f64;
        xstore[BASICLU_DIM] = m as f64;
        xstore[BASICLU_MEMORYL] = 100.0;
        xstore[BASICLU_MEMORYU] = 100.0;
        xstore[BASICLU_MEMORYW] = 100.0;
        (istore, xstore)
    }

    #[test]
    fn load_rejects_bad_hash() {
        let (mut istore, mut xstore) = fresh_store(3);
        istore[0] = 0;
        let (mut li, mut lx) = (vec![0; 4], vec![0.0; 4]);
        let (mut ui, mut ux) = (vec![0; 4], vec![0.0; 4]);
        let (mut wi, mut wx) = (vec![0; 4], vec![0.0; 4]);
        let r = lu_load(
            &mut istore,
            &mut xstore,
            &mut li,
            &mut lx,
            &mut ui,
            &mut ux,
            &mut wi,
            &mut wx,
        );
        assert_eq!(r.err(), Some(BASICLU_ERROR_invalid_store));
    }

    #[test]
    fn load_partitions_istore_and_sets_wend_file_size() {
        let m: LuInt = 3;
        let (mut istore, mut xstore) = fresh_store(m);
        let (mut li, mut lx) = (vec![0; 4], vec![0.0; 4]);
        let (mut ui, mut ux) = (vec![0; 4], vec![0.0; 4]);
        let (mut wi, mut wx) = (vec![0; 4], vec![0.0; 4]);
        {
            let lu = lu_load(
                &mut istore,
                &mut xstore,
                &mut li,
                &mut lx,
                &mut ui,
                &mut ux,
                &mut wi,
                &mut wx,
            )
            .expect("valid store");
            assert_eq!(lu.m, m);
            assert_eq!(lu.colcount_flink.len(), 2 * 3 + 2);
            assert_eq!(lu.wend.len(), 2 * 3 + 1);
            assert_eq!(lu.lbegin_p.len(), 4);
            assert_eq!(lu.iwork0.len(), 3);
            assert_eq!(lu.work0.len(), 3);
            assert_eq!(lu.row_pivot.len(), 3);
        }
        // Wend starts after colcount_flink/blink, rowcount_flink/blink, Wbegin (8+8+8+8+7 entries after istore[1]).
        // nupdate is 0 >= 0 in a zeroed store, so the file size goes to Wend[m].
        assert_eq!(istore[1 + 4 * (2 * 3 + 2) + (2 * 3 + 1) + 3], 100);
        assert_eq!(xstore[BASICLU_DIM] as LuInt, m);
    }

    #[test]
    fn reset_zeroes_work0_only_for_m_bytes() {
        // BASICLU quirk: memset(work0, 0, m*sizeof(lu_int)) clears 4*m bytes.
        // For m = 3 that is one whole double plus the low 4 bytes of the second.
        let m: LuInt = 3;
        let (mut istore, mut xstore) = fresh_store(m);
        let bits = [
            0x1111_1111_2222_2222u64,
            0x3333_3333_4444_4444,
            0x5555_5555_6666_6666,
        ];
        for (k, b) in bits.iter().enumerate() {
            xstore[512 + k] = f64::from_bits(*b);
        }
        let (mut li, mut lx) = (vec![0; 4], vec![0.0; 4]);
        let (mut ui, mut ux) = (vec![0; 4], vec![0.0; 4]);
        let (mut wi, mut wx) = (vec![0; 4], vec![0.0; 4]);
        {
            let mut lu = lu_load(
                &mut istore,
                &mut xstore,
                &mut li,
                &mut lx,
                &mut ui,
                &mut ux,
                &mut wi,
                &mut wx,
            )
            .expect("valid store");
            lu.iwork0.fill(7);
            lu_reset(&mut lu);
            assert_eq!(lu.nupdate, -1);
            assert_eq!(lu.min_colnz, 1);
            assert_eq!(lu.update_cost_denom, 1.0);
            assert!(lu.iwork0.iter().all(|&x| x == 0));
            assert_eq!(lu.wend[2 * 3], 100);
            assert_eq!(lu.work0[0].to_bits(), 0);
            assert_eq!(lu.work0[1].to_bits(), 0x3333_3333_0000_0000);
            assert_eq!(lu.work0[2].to_bits(), bits[2]);
        }
    }

    #[test]
    fn save_writes_scalars_and_status() {
        let m: LuInt = 3;
        let (mut istore, mut xstore) = fresh_store(m);
        let (mut li, mut lx) = (vec![0; 4], vec![0.0; 4]);
        let (mut ui, mut ux) = (vec![0; 4], vec![0.0; 4]);
        let (mut wi, mut wx) = (vec![0; 4], vec![0.0; 4]);
        {
            let mut lu = lu_load(
                &mut istore,
                &mut xstore,
                &mut li,
                &mut lx,
                &mut ui,
                &mut ux,
                &mut wi,
                &mut wx,
            )
            .expect("valid store");
            lu.rank = 2;
            lu.update_cost_numer = 3.0;
            lu.update_cost_denom = 4.0;
            let status = lu_save(&mut lu, BASICLU_WARNING_singular_matrix);
            assert_eq!(status, BASICLU_WARNING_singular_matrix);
        }
        assert_eq!(xstore[BASICLU_RANK], 2.0);
        assert_eq!(
            xstore[BASICLU_STATUS],
            BASICLU_WARNING_singular_matrix as f64
        );
        assert_eq!(xstore[BASICLU_UPDATE_COST], 0.75);
        // user parameters are not written back
        assert_eq!(xstore[BASICLU_MEMORYL], 100.0);
    }
}
