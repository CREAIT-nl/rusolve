//! HiGHS source: ipm/basiclu/lu_update.c
//!
//! Copyright (C) 2016-2018  ERGO-Code
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::lu_dfs::lu_dfs;
use super::lu_file::{lu_file_compress, lu_file_reappend};
use super::lu_garbage_perm::lu_garbage_perm;
use super::lu_internal::Lu;
use super::lu_list::lu_list_swap;
use super::*;

/// lu_update.c: `#define GAP (-1)`
const GAP: LuInt = -1;

/// lu_update.c: `#define FLIP(i) (-(i)-1)`
#[inline]
fn flip(i: LuInt) -> LuInt {
    -i - 1
}

/// HiGHS: basiclu/lu_update.c::find (lines 26-41)
///
/// Find position of index j in index[start..end-1]. If end < 0, then the
/// search stops at the first nonnegative index. Return end if not found.
fn find(j: LuInt, index: &[LuInt], mut start: LuInt, end: LuInt) -> LuInt {
    if end >= 0 {
        while start < end && index[start as usize] != j {
            start += 1;
        }
        start
    } else {
        while index[start as usize] != j && index[start as usize] >= 0 {
            start += 1;
        }
        if index[start as usize] == j {
            start
        } else {
            end
        }
    }
}

/// HiGHS: basiclu/lu_update.c::bfs_path (lines 56-103)
///
/// Find a path from j0 to j0 by a breadth first search. When top < m is
/// returned, then the indices in the path (excluding the final j0) are
/// jlist[top..m-1]. When top == m is returned, then no such path exists.
#[allow(clippy::too_many_arguments)]
fn bfs_path(
    m: LuInt,
    j0: LuInt,
    begin: &[LuInt],
    end: &[LuInt],
    index: &[LuInt],
    jlist: &mut [LuInt],
    marked: &mut [LuInt],
    queue: &mut [LuInt],
) -> LuInt {
    let mut tail: LuInt = 1;
    let mut top = m;
    let mut found = false;
    let mut j = j0;

    queue[0] = j0;
    let mut front: LuInt = 0;
    while front < tail && !found {
        j = queue[front as usize];
        let mut pos = begin[j as usize];
        while pos < end[j as usize] {
            let k = index[pos as usize];
            if k == j0 {
                found = true;
                break;
            }
            if marked[k as usize] >= 0 {
                // parent[k] = j
                marked[k as usize] = flip(j);
                // append to queue
                queue[tail as usize] = k;
                tail += 1;
            }
            pos += 1;
        }
        front += 1;
    }
    if found {
        // build path (j0,..,j)
        while j != j0 {
            top -= 1;
            jlist[top as usize] = j;
            j = flip(marked[j as usize]); // go to parent
            debug_assert!(j >= 0);
        }
        top -= 1;
        jlist[top as usize] = j0;
    }
    for pos in 0..tail {
        marked[queue[pos as usize] as usize] = 0; // reset
    }

    top
}

/// HiGHS: basiclu/lu_update.c::compress_packed (lines 119-172)
///
/// Compress matrix file to reuse memory gaps. Data line 0 <= i < m begins at
/// position begin[i] and ends before the first slot with index[slot] == GAP.
/// begin[m] points to the beginning of unused space at file end.
#[allow(clippy::needless_range_loop)]
fn compress_packed(m: LuInt, begin: &mut [LuInt], index: &mut [LuInt], value: &mut [f64]) -> LuInt {
    let mut nz: LuInt = 0;
    let end = begin[m as usize];

    // Mark the beginning of each nonempty line.
    for i in 0..m as usize {
        let p = begin[i] as usize;
        if index[p] == GAP {
            begin[i] = 0;
        } else {
            debug_assert!(index[p] > GAP);
            begin[i] = index[p]; // temporarily store index here
            index[p] = GAP - i as LuInt - 1; // mark beginning of line i
        }
    }

    // Compress nonempty lines.
    debug_assert!(index[0] == GAP);
    let mut i: LuInt = -1;
    let mut put: LuInt = 1;
    for get in 1..end {
        let g = get as usize;
        if index[g] > GAP {
            // shift entry of line i
            debug_assert!(i >= 0);
            index[put as usize] = index[g];
            value[put as usize] = value[g];
            put += 1;
            nz += 1;
        } else if index[g] < GAP {
            // beginning of line i
            debug_assert!(i == -1);
            i = GAP - index[g] - 1;
            index[put as usize] = begin[i as usize]; // store back
            begin[i as usize] = put;
            value[put as usize] = value[g];
            put += 1;
            nz += 1;
        } else if i >= 0 {
            // line i ended at a gap
            i = -1;
            index[put as usize] = GAP;
            put += 1;
        }
    }
    debug_assert!(i == -1);
    begin[m as usize] = put;
    nz
}

/// HiGHS: basiclu/lu_update.c::permute (lines 192-346)
///
/// Change row-column mappings for columns jlist[0..nswap]. When row i was
/// mapped to column jlist[n], then it will be mapped to column jlist[n+1].
/// When row i was mapped to column jlist[nswap], then it will be mapped to
/// column jlist[0].
///
/// PORT NOTE: C passes `jlist = path + top` with `path = iwork1`. Here the
/// path lives in `this.rowcount_blink`, so the parameter is the offset `top`
/// into it. The entries `jlist[0..=nswap]` are copied out first: `permute`
/// never writes the path, so the copy holds the values C reads.
#[allow(clippy::needless_range_loop)]
fn permute(this: &mut Lu<'_>, jlist_offset: usize, nswap: LuInt) {
    let nswap_u = nswap as usize;
    let jlist: Vec<LuInt> = this.rowcount_blink[jlist_offset..=jlist_offset + nswap_u].to_vec();

    let j0 = jlist[0];
    let jn = jlist[nswap_u];
    let i0 = this.pinv[j0 as usize];
    let in_ = this.pinv[jn as usize];

    debug_assert!(nswap >= 1);
    debug_assert!(this.qinv[i0 as usize] == j0);
    debug_assert!(this.qinv[in_ as usize] == jn);
    debug_assert!(this.row_pivot[i0 as usize] == 0.0);
    debug_assert!(this.col_pivot[j0 as usize] == 0.0);

    // --------------- */
    // Update row file
    // --------------- */
    let begin = this.wbegin[jn as usize]; // keep for later
    let end = this.wend[jn as usize];
    let piv = this.col_pivot[jn as usize];

    let mut n = nswap;
    while n > 0 {
        let j = jlist[n as usize];
        let jprev = jlist[(n - 1) as usize];

        // When row i was indexed by jprev in the row file before,
        // then it is indexed by j now.
        this.wbegin[j as usize] = this.wbegin[jprev as usize];
        this.wend[j as usize] = this.wend[jprev as usize];
        lu_list_swap(this.wflink, this.wblink, j, jprev);

        // That row must have an entry in column j because (jprev,j) is an
        // edge in the augmenting path. This entry becomes a pivot element.
        // If jprev is not the first node in the path, then it has an entry
        // in the row (the old pivot) which becomes an off-diagonal entry now.
        let wpos = find(
            j,
            this.windex,
            this.wbegin[j as usize],
            this.wend[j as usize],
        );
        debug_assert!(wpos < this.wend[j as usize]);
        if n > 1 {
            debug_assert!(jprev != j0);
            this.windex[wpos as usize] = jprev;
            this.col_pivot[j as usize] = this.wvalue[wpos as usize];
            debug_assert!(this.col_pivot[j as usize] != 0.0);
            this.wvalue[wpos as usize] = this.col_pivot[jprev as usize];
        } else {
            debug_assert!(jprev == j0);
            this.col_pivot[j as usize] = this.wvalue[wpos as usize];
            debug_assert!(this.col_pivot[j as usize] != 0.0);
            this.wend[j as usize] -= 1;
            let last = this.wend[j as usize] as usize;
            this.windex[wpos as usize] = this.windex[last];
            this.wvalue[wpos as usize] = this.wvalue[last];
        }
        this.min_pivot = f64::min(this.min_pivot, this.col_pivot[j as usize].abs());
        this.max_pivot = f64::max(this.max_pivot, this.col_pivot[j as usize].abs());
        n -= 1;
    }

    this.wbegin[j0 as usize] = begin;
    this.wend[j0 as usize] = end;
    let wpos = find(
        j0,
        this.windex,
        this.wbegin[j0 as usize],
        this.wend[j0 as usize],
    );
    debug_assert!(wpos < this.wend[j0 as usize]);
    this.windex[wpos as usize] = jn;
    this.col_pivot[j0 as usize] = this.wvalue[wpos as usize];
    debug_assert!(this.col_pivot[j0 as usize] != 0.0);
    this.wvalue[wpos as usize] = piv;
    this.min_pivot = f64::min(this.min_pivot, this.col_pivot[j0 as usize].abs());
    this.max_pivot = f64::max(this.max_pivot, this.col_pivot[j0 as usize].abs());

    // ------------------ */
    // Update column file
    // ------------------ */
    let begin = this.ubegin[i0 as usize]; // keep for later

    for n in 0..nswap_u {
        let i = this.pinv[jlist[n] as usize];
        let inext = this.pinv[jlist[n + 1] as usize];

        // When column j indexed by inext in the column file before,
        // then it is indexed by i now.
        this.ubegin[i as usize] = this.ubegin[inext as usize];

        // That column must have an entry in row i because there is an
        // edge in the augmenting path. This entry becomes a pivot element.
        // There is also an entry in row inext (the old pivot), which now
        // becomes an off-diagonal entry.
        let upos = find(i, this.uindex, this.ubegin[i as usize], -1);
        debug_assert!(upos >= 0);
        this.uindex[upos as usize] = inext;
        this.row_pivot[i as usize] = this.uvalue[upos as usize];
        debug_assert!(this.row_pivot[i as usize] != 0.0);
        this.uvalue[upos as usize] = this.row_pivot[inext as usize];
    }

    this.ubegin[in_ as usize] = begin;
    let upos = find(in_, this.uindex, this.ubegin[in_ as usize], -1);
    debug_assert!(upos >= 0);
    this.row_pivot[in_ as usize] = this.uvalue[upos as usize];
    debug_assert!(this.row_pivot[in_ as usize] != 0.0);
    let mut uend = upos;
    while this.uindex[uend as usize] >= 0 {
        uend += 1;
    }
    this.uindex[upos as usize] = this.uindex[(uend - 1) as usize];
    this.uvalue[upos as usize] = this.uvalue[(uend - 1) as usize];
    this.uindex[(uend - 1) as usize] = -1;

    // -------------------------- */
    // Update row-column mappings
    // -------------------------- */
    let mut n = nswap;
    while n > 0 {
        let j = jlist[n as usize];
        let i = this.pinv[jlist[(n - 1) as usize] as usize];
        this.pinv[j as usize] = i;
        this.qinv[i as usize] = j;
        n -= 1;
    }
    this.pinv[j0 as usize] = in_;
    this.qinv[in_ as usize] = j0;

    for n in 0..=nswap_u {
        let j = jlist[n];
        let i = this.pinv[j as usize];
        debug_assert!(this.row_pivot[i as usize] == this.col_pivot[j as usize]);
    }
}

/// HiGHS: basiclu/lu_update.c::check_consistency (lines 358-411)
///
/// Search for column file entry that is missing or different in row file,
/// and for row file entry that is missing or different in column file.
/// If such an entry is found, then its column and row are returned in *p_col
/// and *p_row. If no such entry exists, then *p_col and *p_row are negative.
///
/// Only called in debug builds (see [`lu_update`]); C compiles it only
/// under `DEBUG_EXTRA`.
#[allow(clippy::needless_range_loop)]
fn check_consistency(this: &Lu<'_>, p_col: &mut LuInt, p_row: &mut LuInt) {
    let m = this.m;

    for i in 0..m as usize {
        // process column file entries
        let j = this.qinv[i];
        let mut pos = this.ubegin[i] as usize;
        loop {
            let ientry = this.uindex[pos];
            if ientry < 0 {
                break;
            }
            let jentry = this.qinv[ientry as usize] as usize;
            let mut wh = this.wbegin[jentry] as usize;
            while (wh as LuInt) < this.wend[jentry] && this.windex[wh] != j {
                wh += 1;
            }
            let found = (wh as LuInt) < this.wend[jentry] && this.wvalue[wh] == this.uvalue[pos];
            if !found {
                *p_col = j;
                *p_row = ientry;
                return;
            }
            pos += 1;
        }
    }
    for j in 0..m as usize {
        // process row file entries
        let i = this.pinv[j];
        let mut pos = this.wbegin[j] as usize;
        while (pos as LuInt) < this.wend[j] {
            let jentry = this.windex[pos];
            let ientry = this.pinv[jentry as usize] as usize;
            let mut wh = this.ubegin[ientry] as usize;
            while this.uindex[wh] >= 0 && this.uindex[wh] != i {
                wh += 1;
            }
            let found = this.uindex[wh] == i && this.uvalue[wh] == this.wvalue[pos];
            if !found {
                *p_col = jentry;
                *p_row = i;
                return;
            }
            pos += 1;
        }
    }
    *p_col = -1;
    *p_row = -1;
}

/// HiGHS: basiclu/lu_update.c::lu_update (lines 430-908)
///
/// Insert spike into U and restore triangularity. If the spiked matrix
/// is permuted triangular, then only permutations are updated. If the
/// spiked matrix is not permuted triangular, then the Forrest-Tomlin
/// update is used and the number of row eta matrices increases by 1.
///
/// Return:
///
///  BASICLU_OK                      update successfully completed
///  BASICLU_REALLOCATE              require more memory in W
///  BASICLU_ERROR_singular_update   new pivot element is zero or < abstol
///
/// PORT NOTE: C `row_reach`/`col_reach` alias `iwork1`/`iwork2`. They are
/// built as owned vectors, because `lu_garbage_perm(this)` needs `&mut Lu`
/// while they are still read afterwards. `lu_garbage_perm` never writes
/// `iwork1`/`iwork2`, so the copies hold the values C reads.
///
/// PORT NOTE: C `pstack = (void *) work1` is passed as `work1` directly,
/// since [`lu_dfs`] takes the `f64` workspace.
///
/// PORT NOTE: `check_consistency` runs under `cfg!(debug_assertions)`. C
/// compiles it only under `DEBUG_EXTRA`; it only asserts, so numerics are
/// unaffected.
#[allow(clippy::too_many_lines)]
pub fn lu_update(this: &mut Lu<'_>, xtbl: f64) -> LuInt {
    let m = this.m;
    let nforrest = this.nforrest;
    let mut unz = this.unz;
    let pad = this.pad;
    let stretch = this.stretch;

    let jpivot = this.btran_for_update;
    let ipivot = this.pinv[jpivot as usize];
    let oldpiv = this.col_pivot[jpivot as usize];

    debug_assert!(nforrest < m);

    // Note: If the singularity test fails or memory is insufficient, then the
    // update is aborted and the user may call this routine a second time.
    // Changes made to data structures in the first call must not violate
    // the logic in the second call.

    // ------- */
    // Prepare
    // ------- */
    // if present, move diagonal element to end of spike
    let mut spike_diag = 0.0_f64;
    let mut have_diag = false;
    let ub_m = this.ubegin[m as usize];
    let mut put = ub_m;
    let mut pos = put;
    loop {
        let i = this.uindex[pos as usize];
        if i < 0 {
            break;
        }
        if i != ipivot {
            this.uindex[put as usize] = i;
            this.uvalue[put as usize] = this.uvalue[pos as usize];
            put += 1;
        } else {
            spike_diag = this.uvalue[pos as usize];
            have_diag = true;
        }
        pos += 1;
    }
    if have_diag {
        this.uindex[put as usize] = ipivot;
        this.uvalue[put as usize] = spike_diag;
    }
    let nz_spike = put - ub_m; // nz excluding diagonal

    let nz_roweta =
        this.rowcount_flink[(nforrest + 1) as usize] - this.rowcount_flink[nforrest as usize];

    // ------------- */
    // Compute pivot
    // ------------- */
    // scatter row eta into work1 and mark positions
    this.marker += 1;
    let mm = this.marker;
    let mut pos = this.rowcount_flink[nforrest as usize];
    while pos < this.rowcount_flink[(nforrest + 1) as usize] {
        let i = this.lindex[pos as usize];
        this.iwork0[i as usize] = mm;
        this.work1[i as usize] = this.lvalue[pos as usize];
        pos += 1;
    }

    // compute newpiv and count intersection
    let mut newpiv = spike_diag;
    let mut intersect: LuInt = 0;
    for pos in ub_m..ub_m + nz_spike {
        let i = this.uindex[pos as usize];
        debug_assert!(i != ipivot);
        if this.iwork0[i as usize] == mm {
            newpiv -= this.uvalue[pos as usize] * this.work1[i as usize];
            intersect += 1;
        }
    }

    // singularity test
    if newpiv == 0.0 || newpiv.abs() < this.abstol {
        return BASICLU_ERROR_singular_update;
    }

    // stability measure
    let piverr = (newpiv - xtbl * oldpiv).abs();

    // ------------ */
    // Insert spike
    // ------------ */
    // calculate bound on file growth
    let mut grow: LuInt = 0;
    for pos in ub_m..ub_m + nz_spike {
        let i = this.uindex[pos as usize];
        debug_assert!(i != ipivot);
        let j = this.qinv[i as usize];
        let jnext = this.wflink[j as usize];
        if this.wend[j as usize] == this.wbegin[jnext as usize] {
            let nz = this.wend[j as usize] - this.wbegin[j as usize];
            grow += nz + 1; // row including spike entry
            // extra room: C `grow += stretch*(nz+1) + pad` (double sum, truncated)
            grow = (grow as f64 + (stretch * ((nz + 1) as f64) + pad as f64)) as LuInt;
        }
    }

    // reallocate if necessary
    let mut room = this.wend[m as usize] - this.wbegin[m as usize];
    if grow > room {
        this.addmem_w = grow - room;
        return BASICLU_REALLOCATE;
    }

    // remove column jpivot from row file
    let mut nz: LuInt = 0;
    let mut pos = this.ubegin[ipivot as usize];
    loop {
        let i = this.uindex[pos as usize];
        if i < 0 {
            break;
        }
        let j = this.qinv[i as usize] as usize;
        let end = this.wend[j];
        this.wend[j] = end - 1;
        let wpos = find(jpivot, this.windex, this.wbegin[j], end);
        debug_assert!(wpos < end);
        this.windex[wpos as usize] = this.windex[(end - 1) as usize];
        this.wvalue[wpos as usize] = this.wvalue[(end - 1) as usize];
        nz += 1;
        pos += 1;
    }
    unz -= nz;

    // erase column jpivot in column file
    let mut pos = this.ubegin[ipivot as usize];
    while this.uindex[pos as usize] >= 0 {
        this.uindex[pos as usize] = GAP;
        pos += 1;
    }

    // set column pointers to spike, chop off diagonal
    this.ubegin[ipivot as usize] = this.ubegin[m as usize];
    this.ubegin[m as usize] += nz_spike;
    let um = this.ubegin[m as usize] as usize;
    this.uindex[um] = GAP;
    this.ubegin[m as usize] += 1;

    // insert spike into row file
    let mut pos = this.ubegin[ipivot as usize];
    loop {
        let i = this.uindex[pos as usize];
        if i < 0 {
            break;
        }
        let j = this.qinv[i as usize] as usize;
        let jnext = this.wflink[j] as usize;
        if this.wend[j] == this.wbegin[jnext] {
            let nz = this.wend[j] - this.wbegin[j];
            // C `room = 1 + stretch*(nz+1) + pad` (double sum, truncated)
            room = (1.0 + stretch * ((nz + 1) as f64) + pad as f64) as LuInt;
            lu_file_reappend(
                j as LuInt,
                m,
                this.wbegin,
                this.wend,
                this.wflink,
                this.wblink,
                this.windex,
                this.wvalue,
                room,
            );
        }
        let end = this.wend[j];
        this.wend[j] += 1;
        this.windex[end as usize] = jpivot;
        this.wvalue[end as usize] = this.uvalue[pos as usize];
        pos += 1;
    }
    unz += nz_spike;

    // insert diagonal
    this.col_pivot[jpivot as usize] = spike_diag;
    this.row_pivot[ipivot as usize] = spike_diag;

    // ------------------ */
    // Test triangularity
    // ------------------ */
    // Some((row_reach, col_reach)) when the spiked matrix stays permuted
    // triangular; None when the Forrest-Tomlin update is needed.
    let triangular: Option<(Vec<LuInt>, Vec<LuInt>)> = if have_diag {
        /*
         * When the spike has a nonzero diagonal element, then the spiked matrix
         * is (symmetrically) permuted triangular if and only if reach(ipivot)
         * does not intersect with the spike pattern except for ipivot. Since
         * reach(ipivot) \ {ipivot} is the structural pattern of the row eta,
         * the matrix is permuted triangular iff the patterns of the row eta
         * and the spike do not intersect.
         */
        let istriangular = intersect == 0;
        if istriangular {
            this.min_pivot = f64::min(this.min_pivot, newpiv.abs());
            this.max_pivot = f64::max(this.max_pivot, newpiv.abs());

            // build row_reach and col_reach in topological order
            let nreach = (nz_roweta + 1) as usize;
            let mut row_reach: Vec<LuInt> = Vec::with_capacity(nreach);
            let mut col_reach: Vec<LuInt> = Vec::with_capacity(nreach);
            row_reach.push(ipivot);
            col_reach.push(jpivot);
            // C: for (n = 1; n < nreach; n++) i = Lindex[pos++]; nreach - 1 == nz_roweta
            let start = this.rowcount_flink[nforrest as usize];
            for pos in start..start + nz_roweta {
                let i = this.lindex[pos as usize];
                row_reach.push(i);
                col_reach.push(this.qinv[i as usize]);
            }
            this.nsymperm_total += 1;
            Some((row_reach, col_reach))
        } else {
            None
        }
    } else {
        /*
         * The spike has a zero diagonal element, so the spiked matrix may only
         * be the *un*symmetric permutation of an upper triangular matrix.
         *
         * Part 1:
         *
         * Find an augmenting path in U[pmap,:] starting from jpivot.
         * An augmenting path is a sequence of column indices such that there
         * is an edge from each node to the next, and an edge from the final
         * node back to jpivot.
         *
         * bfs_path computes such a path in path[top..m-1]. path = iwork1 and
         * reach = iwork2 are the two halves of rowcount_blink.
         */
        let (path, reach) = this.rowcount_blink.split_at_mut(m as usize);
        let top = bfs_path(
            m,
            jpivot,
            this.wbegin,
            this.wend,
            this.windex,
            path,
            this.iwork0,
            reach,
        );
        debug_assert!(top < m - 1);
        debug_assert!(path[top as usize] == jpivot);

        /*
         * Part 2a:
         *
         * For each path index j (except the final one) mark the nodes in
         * reach(j), where the reach is computed in U[pmap,:] without the path
         * edges. If a path index is contained in the reach of an index that
         * comes before it in the path, then U is not permuted triangular.
         *
         * At the same time assemble the combined reach of all path nodes
         * (except the final one) in U[pmap_new,:], where pmap_new is the
         * column-row mapping after applying the permutation associated with
         * the augmenting path. The combined reach is then assembled in
         * topological order in reach[rtop..m-1].
         */
        let mut istriangular = true;
        let mut rtop = m;
        this.marker += 1;
        let mm2 = this.marker;
        let mut t = top;
        while t < m - 1 && istriangular {
            let j = path[t as usize];
            let jnext = path[(t + 1) as usize];
            let wpos = find(
                jnext,
                this.windex,
                this.wbegin[j as usize],
                this.wend[j as usize],
            );
            debug_assert!(wpos < this.wend[j as usize]);
            this.windex[wpos as usize] = j; // take out for a moment
            rtop = lu_dfs(
                j,
                this.wbegin,
                Some(&*this.wend),
                this.windex,
                rtop,
                reach,
                this.work1,
                this.iwork0,
                mm2,
            );
            debug_assert!(reach[rtop as usize] == j);
            reach[rtop as usize] = jnext;
            this.windex[wpos as usize] = jnext; // restore
            istriangular = this.iwork0[jnext as usize] != mm2;
            t += 1;
        }

        /*
         * Part 2b:
         *
         * If the matrix looks triangular so far, then also mark the reach of
         * the final path node, which is reach(jpivot) in U[pmap_new,:].
         * U is then permuted triangular iff the combined reach does not
         * intersect the spike pattern except in the final path index.
         */
        if istriangular {
            let j = path[(m - 1) as usize];
            rtop = lu_dfs(
                j,
                this.wbegin,
                Some(&*this.wend),
                this.windex,
                rtop,
                reach,
                this.work1,
                this.iwork0,
                mm2,
            );
            debug_assert!(reach[rtop as usize] == j);
            reach[rtop as usize] = jpivot;
            this.iwork0[j as usize] -= 1; // unmark for a moment
            let mut pos = this.ubegin[ipivot as usize];
            while this.uindex[pos as usize] >= 0 {
                let i = this.uindex[pos as usize];
                if this.iwork0[this.qinv[i as usize] as usize] == mm2 {
                    istriangular = false;
                }
                pos += 1;
            }
            this.iwork0[j as usize] += 1; // restore
        }

        /*
         * If U is permuted triangular, then permute to zero-free diagonal.
         * Set up row_reach[0..nreach-1] and col_reach[0..nreach-1] for
         * updating the permutations below. The column reach is the combined
         * reach of the path nodes, copied out before the permutation. The row
         * reach is given through pmap, which permute changes, so it is read
         * after the permutation.
         */
        if istriangular {
            let nswap = m - top - 1;
            debug_assert!(reach[rtop as usize] == jpivot);
            let cols: Vec<LuInt> = reach[rtop as usize..m as usize].to_vec();
            permute(this, top as usize, nswap);
            unz -= 1;
            let rows: Vec<LuInt> = cols.iter().map(|&c| this.pinv[c as usize]).collect();
            Some((rows, cols))
        } else {
            None
        }
    };

    // --------------------- */
    // Forrest-Tomlin update
    // --------------------- */
    let (row_reach, col_reach) = match triangular {
        Some(rc) => rc,
        None => {
            // remove row ipivot from column file
            let mut pos = this.wbegin[jpivot as usize];
            while pos < this.wend[jpivot as usize] {
                let j = this.windex[pos as usize];
                debug_assert!(j != jpivot);
                let mut wp: LuInt = -1;
                let mut end = this.ubegin[this.pinv[j as usize] as usize];
                loop {
                    let i = this.uindex[end as usize];
                    if i < 0 {
                        break;
                    }
                    if i == ipivot {
                        wp = end;
                    }
                    end += 1;
                }
                debug_assert!(wp >= 0);
                this.uindex[wp as usize] = this.uindex[(end - 1) as usize];
                this.uvalue[wp as usize] = this.uvalue[(end - 1) as usize];
                this.uindex[(end - 1) as usize] = -1;
                unz -= 1;
                pos += 1;
            }

            // remove row ipivot from row file
            this.wend[jpivot as usize] = this.wbegin[jpivot as usize];

            // replace pivot
            this.col_pivot[jpivot as usize] = newpiv;
            this.row_pivot[ipivot as usize] = newpiv;
            this.min_pivot = f64::min(this.min_pivot, newpiv.abs());
            this.max_pivot = f64::max(this.max_pivot, newpiv.abs());

            // drop zeros from row eta; update max entry of row etas
            let mut nz: LuInt = 0;
            let mut put = this.rowcount_flink[nforrest as usize];
            let mut max_eta = 0.0_f64;
            let mut pos = put;
            while pos < this.rowcount_flink[(nforrest + 1) as usize] {
                if this.lvalue[pos as usize] != 0.0 {
                    max_eta = f64::max(max_eta, this.lvalue[pos as usize].abs());
                    this.lindex[put as usize] = this.lindex[pos as usize];
                    this.lvalue[put as usize] = this.lvalue[pos as usize];
                    put += 1;
                    nz += 1;
                }
                pos += 1;
            }
            this.rowcount_flink[(nforrest + 1) as usize] = put;
            this.rnz += nz;
            this.max_eta = f64::max(this.max_eta, max_eta);

            // prepare permutation update
            this.nforrest += 1;
            this.nforrest_total += 1;
            (vec![ipivot], vec![jpivot])
        }
    };

    // ------------------- */
    // Update permutations
    // ------------------- */
    let nreach = row_reach.len() as LuInt;
    if this.pivotlen + nreach > 2 * m {
        lu_garbage_perm(this);
    }

    // append row indices row_reach[0..nreach-1] to end of pivot sequence
    let base = this.pivotlen as usize;
    for (n, &r) in row_reach.iter().enumerate() {
        this.colcount_blink[base + n] = r; // pivotrow
    }

    // append col indices col_reach[0..nreach-1] to end of pivot sequence
    for (n, &c) in col_reach.iter().enumerate() {
        this.colcount_flink[base + n] = c; // pivotcol
    }

    this.pivotlen += nreach;

    // -------- */
    // Clean up
    // -------- */
    // compress U if used memory is shrinked sufficiently
    let used = this.ubegin[m as usize];
    if ((used - unz - m) as f64) > this.compress_thres * (used as f64) {
        let nz = compress_packed(m, this.ubegin, this.uindex, this.uvalue);
        debug_assert_eq!(nz, unz);
    }

    // compress W if used memory is shrinked sufficiently
    let used = this.wbegin[m as usize];
    // C `need = Unz + stretch*Unz + m*pad` (double sum, truncated)
    let need = (unz as f64 + stretch * (unz as f64) + (m * pad) as f64) as LuInt;
    if ((used - need) as f64) > this.compress_thres * (used as f64) {
        let nz = lu_file_compress(
            m,
            this.wbegin,
            this.wend,
            this.wflink,
            this.windex,
            this.wvalue,
            stretch,
            pad,
        );
        debug_assert_eq!(nz, unz);
    }

    this.pivot_error = piverr / (1.0 + newpiv.abs());
    this.unz = unz;
    this.btran_for_update = -1;
    this.ftran_for_update = -1;
    this.update_cost_numer += nz_roweta as f64;
    this.nupdate += 1;
    this.nupdate_total += 1;

    if cfg!(debug_assertions) {
        let mut col: LuInt = 0;
        let mut row: LuInt = 0;
        check_consistency(this, &mut col, &mut row);
        debug_assert!(col < 0 && row < 0);
    }

    BASICLU_OK
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_respects_bounded_and_sentinel_scans() {
        // end >= 0: bounded scan, returns end when not found
        assert_eq!(find(4, &[3, 1, 4, 2], 0, 2), 2);
        assert_eq!(find(4, &[3, 1, 4, 2], 0, 4), 2);
        assert_eq!(find(2, &[3, 1, 4, 2], 0, 4), 3);
        // end < 0: scan stops at the first negative (gap) entry, returns end
        assert_eq!(find(1, &[3, 1, -1, 4], 0, -1), 1);
        assert_eq!(find(4, &[3, 1, -1, 4], 0, -1), -1);
    }

    #[test]
    fn bfs_finds_cycle_path_and_resets_marks() {
        // Graph on 3 nodes with the cycle 0 -> 1 -> 2 -> 0.
        let begin = [0, 1, 2];
        let end = [1, 2, 3];
        let index = [1, 2, 0];
        let mut jlist = [0; 3];
        let mut marked = [0; 3];
        let mut queue = [0; 3];
        let top = bfs_path(
            3,
            0,
            &begin,
            &end,
            &index,
            &mut jlist,
            &mut marked,
            &mut queue,
        );
        // path (0, 1, 2): jlist[top..m] = [0, 1, 2]
        assert_eq!(top, 0);
        assert_eq!(jlist, [0, 1, 2]);
        // marks are reset to zero after the search
        assert_eq!(marked, [0, 0, 0]);
    }

    #[test]
    fn compress_packed_moves_lines_into_gap_free_prefix() {
        // Line 0 = {5, 7} at slots 1..2, gap at 3; line 1 = {9} at slot 4,
        // gap at 5; unused space from slot 6.
        let mut begin = [1, 4, 6];
        let mut index = [-1, 5, 7, -1, 9, -1];
        let mut value = [0.0, 50.0, 70.0, 0.0, 90.0, 0.0];
        let nz = compress_packed(2, &mut begin, &mut index, &mut value);
        assert_eq!(nz, 3);
        assert_eq!(begin, [1, 4, 6]);
        assert_eq!(&index[1..], &[5, 7, -1, 9, -1]);
        assert_eq!(value[1], 50.0);
        assert_eq!(value[2], 70.0);
        assert_eq!(value[4], 90.0);
    }
}
