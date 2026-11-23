//! HiGHS source: ipm/ipx/maxvolume.h, maxvolume.cc
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.

use super::basis::{BASIC, Basis, NONBASIC};
use super::control::{Control, sci2};
use super::indexed_vector::{IndexedVector, for_each_nonzero, for_each_nonzero_mut};
use super::sparse_matrix::dot_column;
use super::timer::Timer;
use super::utils::sortperm;
use super::*;
use crate::util::highs_utils::{cpp_max, cpp_min};
use std::rc::Rc;

/// HiGHS: ipx/maxvolume.cc::Maxvolume::Slice (lines 94-104)
struct Slice {
    colscale: Vector,
    invscale_basic: Vector,
    tblrow_used: Vec<bool>,
    colweights: Vector,
    lhs: IndexedVector,
    row: IndexedVector,
    work: Vector,
}

impl Slice {
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::Slice::Slice (lines 95-97)
    fn new(m: Int, n: Int) -> Self {
        Slice {
            colscale: Vector::new(n + m),
            invscale_basic: Vector::new(m),
            tblrow_used: vec![false; m as usize],
            colweights: Vector::new(n + m),
            lhs: IndexedVector::new(m),
            row: IndexedVector::new(n + m),
            work: Vector::new(m),
        }
    }
}

/// HiGHS: ipx/maxvolume.h::Maxvolume
pub struct Maxvolume {
    control_: Rc<Control>,
    updates_: Int,
    skipped_: Int,
    passes_: Int,
    slices_: Int,
    volinc_: f64,
    time_: f64,
    tblnnz_: Int,
    tblmax_: f64,
    frobnorm_squared_: f64,
}

impl Maxvolume {
    pub const kPivotZeroTol: f64 = 1e-7;

    /// HiGHS: ipx/maxvolume.cc::Maxvolume::Maxvolume (line 10)
    pub fn new(control: Rc<Control>) -> Self {
        Maxvolume {
            control_: control,
            updates_: 0,
            skipped_: 0,
            passes_: 0,
            slices_: 0,
            volinc_: 0.0,
            time_: 0.0,
            tblnnz_: 0,
            tblmax_: 0.0,
            frobnorm_squared_: 0.0,
        }
    }

    /// HiGHS: ipx/maxvolume.cc::Maxvolume::RunSequential (lines 12-92)
    ///
    /// IPX quirk: the C++ `colscale` pointer may be null (`colscale ? ... : 1.0`);
    /// every call site passes a non-null array, so the Rust slice is always present.
    pub fn run_sequential(&mut self, colscale: &[f64], basis: &mut Basis) -> Int {
        let model = basis.model_rc();
        let m = model.rows();
        let n = model.cols();
        let mut ftran = IndexedVector::new(m);
        let timer = Timer::new(0.0);
        let mut errflag: Int = 0;

        let maxpasses = self.control_.maxpasses();
        let volumetol = cpp_max(self.control_.volume_tol(), 1.0);

        // Maintain a copy of the inverse scaling factors of basic variables.
        // If a variable has status BASIC_FREE, its inverse scaling factor is
        // assumed zero, regardless of the entry in colscale. This ensures that
        // BASIC_FREE variables are never pivoted out of the basis.
        let mut invscale_basic = Vector::new(m);
        for p in 0..m {
            let j = basis.at(p);
            if basis.status_of(j) == BASIC {
                invscale_basic[p as usize] = 1.0 / colscale[j as usize];
            }
        }

        self.reset();
        while self.passes_ < maxpasses || maxpasses < 0 {
            self.tblnnz_ = 0;
            self.tblmax_ = 0.0;
            self.frobnorm_squared_ = 0.0;
            let mut updates_last: Int = 0; // # basis updates in this pass
            let mut candidates = sortperm(n + m, Some(colscale), false);
            while let Some(&j) = candidates.last() {
                let ju = j as usize;
                let dj = colscale[ju];
                if dj == 0.0 {
                    // all remaining columns have scaling factor 0
                    break;
                }
                if basis.status_of(j) != NONBASIC {
                    candidates.pop();
                    continue;
                }
                errflag = self.control_.interrupt_check(-1);
                if errflag != 0 {
                    break;
                }
                basis.solve_for_update_lhs(j, &mut ftran);
                let mut pmax: Int = -1;
                let mut vmax: f64 = 0.0;
                for_each_nonzero(&ftran, |p, x| {
                    let v = x.abs() * invscale_basic[p as usize] * dj;
                    if v > vmax {
                        vmax = v;
                        pmax = p;
                    }
                    if v != 0.0 {
                        self.tblnnz_ += 1;
                    }
                    self.frobnorm_squared_ += v * v;
                });
                self.tblmax_ = cpp_max(self.tblmax_, vmax);
                if vmax <= volumetol {
                    self.skipped_ += 1;
                    candidates.pop();
                    continue;
                }

                let jb = basis.at(pmax);
                debug_assert!(basis.status_of(jb) == BASIC);
                let mut exchanged = false;
                errflag = basis.exchange_if_stable(jb, j, ftran.get(pmax), -1, &mut exchanged);
                if errflag != 0 {
                    break;
                }
                if !exchanged {
                    // factorization was unstable, try again
                    continue;
                }
                invscale_basic[pmax as usize] = 1.0 / dj;
                updates_last += 1;
                self.volinc_ += vmax.log2();
                candidates.pop();
            }
            self.updates_ += updates_last;
            self.passes_ += 1;
            if updates_last == 0 || errflag != 0 {
                break;
            }
        }
        self.time_ = timer.elapsed();
        errflag
    }

    /// HiGHS: ipx/maxvolume.cc::Maxvolume::RunHeuristic (lines 106-151)
    pub fn run_heuristic(&mut self, colscale: &[f64], basis: &mut Basis) -> Int {
        let model = basis.model_rc();
        let m = model.rows();
        let n = model.cols();
        let mut slice = Slice::new(m, n);
        let mut errflag: Int = 0;
        let timer = Timer::new(0.0);

        self.reset();
        // IPX quirk: C++ `5 + std::max((long)(...), 0l)` is computed in 64-bit
        // `long` and then narrowed to `Int` on assignment; `as Int` wraps the same way.
        let rows_per_slice = self.control_.rows_per_slice();
        let mut num_slices: Int = (5 + cpp_max(i64::from(m / rows_per_slice), 0)) as Int;
        num_slices = cpp_min(num_slices, m);

        // Maintain a copy of the inverse scaling factors of basic variables.
        for p in 0..m {
            let j = basis.at(p);
            if basis.status_of(j) == BASIC {
                slice.invscale_basic[p as usize] = 1.0 / colscale[j as usize];
                debug_assert!(slice.invscale_basic[p as usize].is_finite());
            }
        }

        // We need to work with a copy of the column scaling factors, in which
        // skipped columns are set to zero (and not scanned again).
        for j in 0..n + m {
            if basis.status_of(j) == NONBASIC {
                slice.colscale[j as usize] = colscale[j as usize];
            }
        }

        // Split tableau matrix into num_slices row slices. In each call to
        // Driver() there is exactly one row from each slice that has
        // tblrow_used[p] == 1, while the other rows in the slice have
        // tblrow_used[p] == 0.
        let perm = sortperm(m, Some(&slice.invscale_basic[..]), false);
        for s in 0..num_slices {
            for i in 0..m {
                slice.tblrow_used[perm[i as usize] as usize] = i % num_slices == s;
            }
            errflag = self.driver(basis, &mut slice);
            if errflag != 0 {
                break;
            }
        }

        self.time_ = timer.elapsed();
        self.passes_ = -1;
        self.slices_ = num_slices;
        errflag
    }

    /// HiGHS: ipx/maxvolume.cc::Maxvolume::updates (line 153)
    pub fn updates(&self) -> Int {
        self.updates_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::skipped (line 154)
    pub fn skipped(&self) -> Int {
        self.skipped_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::passes (line 155)
    pub fn passes(&self) -> Int {
        self.passes_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::slices (line 156)
    pub fn slices(&self) -> Int {
        self.slices_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::volinc (line 157)
    pub fn volinc(&self) -> f64 {
        self.volinc_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::time (line 158)
    pub fn time(&self) -> f64 {
        self.time_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::tblnnz (line 159)
    pub fn tblnnz(&self) -> Int {
        self.tblnnz_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::tblmax (line 160)
    pub fn tblmax(&self) -> f64 {
        self.tblmax_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::frobnorm_squared (line 161)
    pub fn frobnorm_squared(&self) -> f64 {
        self.frobnorm_squared_
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::Reset (lines 163-173)
    fn reset(&mut self) {
        self.updates_ = 0;
        self.skipped_ = 0;
        self.passes_ = 0;
        self.slices_ = 0;
        self.volinc_ = 0.0;
        self.time_ = 0.0;
        self.tblnnz_ = 0;
        self.tblmax_ = 0.0;
        self.frobnorm_squared_ = 0.0;
    }
    /// HiGHS: ipx/maxvolume.cc::Maxvolume::Driver (lines 200-318)
    fn driver(&mut self, basis: &mut Basis, slice: &mut Slice) -> Int {
        let model = basis.model_rc();
        let m = model.rows();
        let n = model.cols();
        let ai = model.AI();
        let mut errflag: Int = 0;

        let volumetol = cpp_max(self.control_.volume_tol(), 1.0);
        let maxskip = self.control_.maxskip_updates();

        let Slice {
            colscale,
            invscale_basic,
            tblrow_used,
            colweights,
            lhs,
            row,
            work,
        } = slice;

        // Compute column weights.
        for p in 0..m {
            let pu = p as usize;
            work[pu] = if tblrow_used[pu] {
                invscale_basic[pu]
            } else {
                0.0
            };
        }
        // IPX quirk: C++ aliases rhs and lhs in SolveDense(work, work, 'T');
        // BASICLU copies rhs before writing lhs, so cloning rhs is exact.
        let rhs = work.clone();
        basis.solve_dense(&rhs, work, b'T');
        for j in 0..n + m {
            let ju = j as usize;
            if colscale[ju] != 0.0 {
                debug_assert!(basis.status_of(j) == NONBASIC);
                let sum = dot_column(ai, j, &work[..]);
                colweights[ju] = sum * colscale[ju];
            } else {
                colweights[ju] = 0.0;
            }
        }

        let mut candidates: Vec<Int> = Vec::new();
        let mut skipped: Int = 0;
        loop {
            // Pick column with maximum weight.
            if candidates.is_empty() {
                candidates = find_largest(colweights);
                debug_assert!(!candidates.is_empty());
            }
            let jn = candidates[candidates.len() - 1];
            debug_assert!(jn >= 0 && jn < n + m);
            let weight = colweights[jn as usize];
            if weight == 0.0 {
                break;
            }
            debug_assert!(basis.status_of(jn) == NONBASIC);
            debug_assert!(colscale[jn as usize] > 0.0);

            errflag = self.control_.interrupt_check(-1);
            if errflag != 0 {
                break;
            }

            // Find maximum scaled FTRAN entry.
            basis.solve_for_update_lhs(jn, lhs);
            let pmax = self.scale_ftran(colscale[jn as usize], invscale_basic, lhs);
            let scaled_pivot = lhs.get(pmax);
            let vmax = scaled_pivot.abs();

            // Skip the column if exchange does not increase volume enough.
            if vmax <= volumetol {
                colweights[jn as usize] = 0.0;
                colscale[jn as usize] = 0.0;
                candidates.pop();
                skipped += 1;
                if skipped > maxskip && maxskip >= 0 {
                    break;
                }
                continue;
            }

            // Recompute column weight from FTRAN.
            let mut weight_recomp: f64 = 0.0;
            for_each_nonzero(lhs, |p, x| {
                if tblrow_used[p as usize] {
                    weight_recomp += x;
                }
            });
            debug_assert!(weight_recomp.is_finite());

            // Update basis.
            let jb = basis.at(pmax);
            basis.tableau_row(jb, lhs, row, true);
            let pivot = row.get(jn);
            if pivot.abs() < 1e-3 && self.control_.debug_enabled(3) {
                self.control_
                    .debug(3, &format!(" |pivot| {}(maxvolume)\n", sci2(pivot.abs())));
            }
            let mut exchanged = false;
            errflag = basis.exchange_if_stable(jb, jn, pivot, 0, &mut exchanged);
            if errflag != 0 {
                break;
            }
            if !exchanged {
                // factorization was unstable, try again
                continue;
            }
            self.updates_ += 1;
            self.volinc_ += vmax.log2();

            // Update colscale and invscale_basic.
            let dn = colscale[jn as usize];
            let dbinv = invscale_basic[pmax as usize];
            debug_assert!(colscale[jb as usize] == 0.0);
            colscale[jb as usize] = 1.0 / invscale_basic[pmax as usize];
            invscale_basic[pmax as usize] = 1.0 / colscale[jn as usize];
            colscale[jn as usize] = 0.0;
            debug_assert!(colscale[jb as usize].is_finite());
            debug_assert!(invscale_basic[pmax as usize].is_finite());

            // Update column weights. C++ `bool` promotes to 0.0 / 1.0 here.
            let used: f64 = if tblrow_used[pmax as usize] { 1.0 } else { 0.0 };
            let alpha = (used - weight_recomp) / (dn * pivot);
            debug_assert!(alpha.is_finite());
            for_each_nonzero(row, |j, x| {
                colweights[j as usize] += alpha * x * colscale[j as usize];
            });
            colweights[jb as usize] = used + alpha / dbinv;
            colweights[jn as usize] = 0.0;
            candidates.clear();
        }

        self.skipped_ += skipped;
        errflag
    }

    /// HiGHS: ipx/maxvolume.cc::Maxvolume::ScaleFtran (lines 320-335)
    fn scale_ftran(
        &mut self,
        colscale_jn: f64,
        invscale_basic: &Vector,
        ftran: &mut IndexedVector,
    ) -> Int {
        let mut vmax: f64 = 0.0;
        let mut pmax: Int = 0;
        for_each_nonzero_mut(ftran, |p, pivot| {
            let scaled_pivot = *pivot * colscale_jn * invscale_basic[p as usize];
            let v = scaled_pivot.abs();
            if v > vmax && pivot.abs() > Self::kPivotZeroTol {
                vmax = v;
                pmax = p;
            }
            *pivot = scaled_pivot;
        });
        pmax
    }
}

/// HiGHS: ipx/maxvolume.cc::FindLargest (static, lines 177-198)
///
/// Returns the indices of the largest and second largest entry in `weights`.
/// The index of the largest entry is at the back. Ties go to the first index
/// (strict `>` comparisons), as in the C++ source.
fn find_largest(weights: &Vector) -> Vec<Int> {
    let n = weights.size();
    let mut jmax: Int = 0; // index of largest element
    let mut jmax2: Int = 0; // index of second largest element
    let mut wmax: f64 = 0.0;
    let mut wmax2: f64 = 0.0;

    for j in 0..n {
        let w = weights[j as usize].abs();
        if w > wmax {
            wmax2 = wmax;
            wmax = w;
            jmax2 = jmax;
            jmax = j;
        } else if w > wmax2 {
            wmax2 = w;
            jmax2 = j;
        }
    }
    vec![jmax2, jmax]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_largest_picks_two_largest_with_first_index_on_ties() {
        // |w| = [1, 5, 3, 5, 0]. jmax=1 (first 5), jmax2=3 (second 5 wins the
        // second slot because 5 > 3 and it is not > wmax=5).
        let w = Vector::from_slice(&[1.0, -5.0, 3.0, 5.0, 0.0]);
        assert_eq!(find_largest(&w), vec![3, 1]);
    }
}
