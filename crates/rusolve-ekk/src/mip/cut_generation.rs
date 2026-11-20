//! HiGHS source: mip/HighsCutGeneration.h, mip/HighsCutGeneration.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.
//!
//! PORT NOTE: C++ `assert(..)` calls are not mirrored (release / NDEBUG semantics). The one
//! `#ifndef NDEBUG` check of `cmirCutGenerationHeuristic` is mirrored with `cfg!(debug_assertions)`.
//! The `#if 0` blocks of the C++ file are not mirrored.

use crate::prelude::*;
use crate::util::highs_hash::{frexp, hash_u32_pair};
use crate::util::highs_integers::HighsIntegers;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// HiGHS: HighsCutGeneration.h::HighsCutGeneration (lines 31-114)
///
/// Removed C++ members:
/// * C++: `const HighsLpRelaxation & lpRelaxation`: removed; methods that use it take
///   `lp_relaxation: &HighsLpRelaxation` (detached by the caller, the separator's `&mut` LP).
///   `lpRelaxation.getMipSolver()` is the `ms` parameter.
/// * C++: `HighsCutPool & cutpool`: kept as the index [`HighsCutGeneration::cutpool`] into
///   `ms.mipdata_.cutpools`; `cutpool.addCut(..)` is
///   `detach!(ms.mipdata_.cutpools[self.cutpool], |p| p.add_cut(ms, ..))`.
/// * C++: `double * vals` / `HighsInt * inds`: removed. They always point at the data of one
///   of the vectors the public entry points receive (`inds_` / `vals_`, `proofinds` /
///   `proofvals`, `indsCheck_` / `valsCheck_` under debug, or `tmpInds` / `tmpVals` in
///   `tryGenerateCut`). The public entry points take only those vectors and pass
///   `&mut v[..]` to the private helpers, which take `vals: &mut [f64]` (and
///   `inds: &mut [HighsInt]` where the C++ body reads `inds`). Every C++ re-point
///   `inds = X.data(); vals = Y.data();` is "the slices passed from here on are `X` / `Y`".
#[derive(Clone, Default, Debug)]
pub struct HighsCutGeneration {
    /// rusolve (mip_tableau_loosen): no base-row length cap
    pub loosen_len: bool,
    /// C++: `HighsRandom randgen`
    pub randgen: HighsRandom,
    /// C++: `std::vector<HighsInt> cover`
    pub cover: Vec<HighsInt>,
    /// C++: `HighsCDouble coverweight`
    pub coverweight: HighsCDouble,
    /// C++: `HighsCDouble lambda`
    pub lambda: HighsCDouble,
    /// C++: `std::vector<double> upper`
    pub upper: Vec<f64>,
    /// C++: `std::vector<double> solval`
    pub solval: Vec<f64>,
    /// C++: `std::vector<uint8_t> complementation`
    pub complementation: Vec<u8>,
    /// C++: `std::vector<uint8_t> isintegral`
    pub isintegral: Vec<u8>,
    /// C++: `const double feastol`
    pub feastol: f64,
    /// C++: `const double epsilon`
    pub epsilon: f64,
    /// C++: `HighsCDouble rhs`
    pub rhs: HighsCDouble,
    /// C++: `bool integralSupport`
    pub integral_support: bool,
    /// C++: `bool integralCoefficients`
    pub integral_coefficients: bool,
    /// C++: `HighsInt rowlen`
    pub rowlen: HighsInt,
    /// C++: `double initialScale`
    pub initial_scale: f64,
    /// C++: `std::vector<HighsInt> integerinds`
    pub integerinds: Vec<HighsInt>,
    /// C++: `std::vector<double> deltas`
    pub deltas: Vec<f64>,
    /// C++: `std::vector<double> tmpVals`
    pub tmp_vals: Vec<f64>,
    /// C++: `std::vector<HighsInt> tmpInds`
    pub tmp_inds: Vec<HighsInt>,
    /// C++: `std::vector<uint8_t> tmpComplementation`
    pub tmp_complementation: Vec<u8>,
    /// C++: `std::vector<double> tmpSolval`
    pub tmp_solval: Vec<f64>,
    /// C++: `HighsCutPool & cutpool` (index into `ms.mipdata_.cutpools`; see the struct doc).
    pub cutpool: usize,
    /// Scratch for `cmir_cut_generation_heuristic`: `(vals[j], solval[j], upper[j])` for
    /// each `j` in `integerinds`, gathered once per call (none of the three changes while
    /// the scale loops run). Not in the C++; it stands in for the indexed reads.
    pub cmir_gather: Vec<(f64, f64, f64)>,
}

impl HighsCutGeneration {
    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::HighsCutGeneration (lines 16-23)
    /// C++: `HighsCutGeneration(const HighsLpRelaxation & lpRelaxation, HighsCutPool & cutpool)`
    // PORT NOTE: `cutpool` is an index into `ms.mipdata_.cutpools` (attached here). The C++
    // seed is the `int64_t` sum `random_seed + getNumLpIterations() + getNumCuts()` converted
    // to `HighsUInt` (mod 2^32), which `as HighsUInt` reproduces.
    pub fn new(ms: &HighsMipSolver, lp_relaxation: &HighsLpRelaxation, cutpool: usize) -> Self {
        let seed = ms.options_mip_.random_seed as i64
            + lp_relaxation.get_num_lp_iterations()
            + ms.mipdata_.cutpools[cutpool].get_num_cuts() as i64;
        HighsCutGeneration {
            randgen: HighsRandom::new(seed as HighsUInt),
            feastol: ms.mipdata_.feastol,
            epsilon: ms.mipdata_.epsilon,
            cutpool,
            ..Default::default()
        }
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::determineCover (lines 25-129)
    /// C++: `bool determineCover(bool lpSol = true)`
    /// C++ default arguments: `lp_sol = true`
    // PORT NOTE: `lpRelaxation.getMipSolver().mipdata_->nodequeue` is `ms.mipdata_.nodequeue`
    // (attached).
    pub fn determine_cover(
        &mut self,
        ms: &HighsMipSolver,
        vals: &mut [f64],
        inds: &mut [HighsInt],
        lp_sol: bool,
    ) -> bool {
        if self.rhs <= 10.0 * self.feastol {
            return false;
        }

        self.cover.clear();

        for j in 0..self.rowlen {
            let ju = j as usize;
            if self.isintegral[ju] == 0 {
                continue;
            }

            if lp_sol && self.solval[ju] <= self.feastol {
                continue;
            }

            self.cover.push(j);
        }

        let max_cover_size = self.cover.len();
        let mut coversize: usize = 0;
        let r: HighsInt = self.randgen.integer();
        self.coverweight = HighsCDouble::from(0.0);
        if lp_sol {
            let feastol = self.feastol;
            // take all variables that sit at their upper bound always into the cover
            {
                let solval = &self.solval;
                let upper = &self.upper;
                coversize = cpp_partition(&mut self.cover[..], |&j| {
                    solval[j as usize] >= upper[j as usize] - feastol
                });
            }

            for i in 0..coversize {
                let j = self.cover[i] as usize;
                self.coverweight += vals[j] * self.upper[j];
            }

            // sort the remaining variables by the contribution to the rows activity in
            // the current solution
            let solval = &self.solval;
            let upper = &self.upper;
            let vals_r: &[f64] = &*vals;
            let inds_r: &[HighsInt] = &*inds;
            pdqsort_by(&mut self.cover[coversize..max_cover_size], |&i, &j| {
                let (iu, ju) = (i as usize, j as usize);
                if upper[iu] < 1.5 && upper[ju] > 1.5 {
                    return true;
                }
                if upper[iu] > 1.5 && upper[ju] < 1.5 {
                    return false;
                }

                let contribution_a = solval[iu] * vals_r[iu];
                let contribution_b = solval[ju] * vals_r[ju];

                if contribution_a > contribution_b + feastol {
                    return true;
                }
                if contribution_a < contribution_b - feastol {
                    return false;
                }
                // for equal contributions take the larger coefficients first
                // because this makes some of the lifting functions more likely
                // to generate a facet
                // if the value is equal too, choose a random tiebreaker based
                // on hashing the column index and the current number of pool
                // cuts
                if (vals_r[iu] - vals_r[ju]).abs() <= feastol {
                    return hash_u32_pair(inds_r[iu] as u32, r as u32)
                        > hash_u32_pair(inds_r[ju] as u32, r as u32);
                }
                vals_r[iu] > vals_r[ju]
            });
        } else {
            // the current solution
            let nodequeue = &ms.mipdata_.nodequeue;
            let feastol = self.feastol;
            let solval = &self.solval;
            let complementation = &self.complementation;
            let inds_r: &[HighsInt] = &*inds;

            pdqsort_by(&mut self.cover[coversize..max_cover_size], |&i, &j| {
                let (iu, ju) = (i as usize, j as usize);
                if solval[iu] > feastol && solval[ju] <= feastol {
                    return true;
                }
                if solval[iu] <= feastol && solval[ju] > feastol {
                    return false;
                }

                let num_nodes_a: i64 = if complementation[iu] != 0 {
                    nodequeue.num_nodes_down(inds_r[iu])
                } else {
                    nodequeue.num_nodes_up(inds_r[iu])
                };

                let num_nodes_b: i64 = if complementation[ju] != 0 {
                    nodequeue.num_nodes_down(inds_r[ju])
                } else {
                    nodequeue.num_nodes_up(inds_r[ju])
                };

                if num_nodes_a > num_nodes_b {
                    return true;
                }
                if num_nodes_a < num_nodes_b {
                    return false;
                }

                hash_u32_pair(inds_r[iu] as u32, r as u32)
                    > hash_u32_pair(inds_r[ju] as u32, r as u32)
            });
        }

        let minlambda = cpp_max(
            10.0 * self.feastol,
            self.feastol * f64::from(self.rhs).abs(),
        );

        while coversize != max_cover_size {
            let lambda_local = f64::from(self.coverweight - self.rhs);
            if lambda_local > minlambda {
                break;
            }

            let j = self.cover[coversize] as usize;
            self.coverweight += vals[j] * self.upper[j];
            coversize += 1;
        }
        if coversize == 0 {
            return false;
        }

        self.coverweight.renormalize();
        self.lambda = self.coverweight - self.rhs;

        if self.lambda <= minlambda {
            return false;
        }

        self.cover.truncate(coversize);

        true
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::separateLiftedKnapsackCover (lines 131-214)
    /// C++: `void separateLiftedKnapsackCover()`
    #[allow(clippy::needless_range_loop)]
    pub fn separate_lifted_knapsack_cover(&mut self, vals: &mut [f64]) {
        let coversize = self.cover.len();

        let mut s = vec![0.0f64; coversize];
        let mut coverflag = vec![0i8; self.rowlen as usize];
        pdqsort_branchless_by(&mut self.cover[..], |&a, &b| {
            vals[a as usize] > vals[b as usize]
        });

        let mut abartmp = HighsCDouble::from(vals[self.cover[0] as usize]);
        let mut sigma = self.lambda;
        for i in 1..coversize {
            let vc = vals[self.cover[i] as usize];
            let delta = abartmp - vc;
            let kdelta = (i as f64) * delta;
            if f64::from(kdelta) < f64::from(sigma) {
                abartmp = HighsCDouble::from(vc);
                sigma -= kdelta;
            } else {
                abartmp -= sigma * (1.0 / (i as f64));
                sigma = HighsCDouble::from(0.0);
                break;
            }
        }

        if f64::from(sigma) > 0.0 {
            abartmp = self.rhs / (coversize as f64);
        }

        let abar = f64::from(abartmp);

        let mut sum = HighsCDouble::from(0.0);
        let mut cplussize: HighsInt = 0;
        for i in 0..coversize {
            let ci = self.cover[i] as usize;
            sum += cpp_min(abar, vals[ci]);
            s[i] = f64::from(sum);

            if vals[ci] > abar + self.feastol {
                cplussize += 1;
                coverflag[ci] = 1;
            } else {
                coverflag[ci] = -1;
            }
        }
        let mut halfintegral = false;

        let feastol = self.feastol;
        let epsilon = self.epsilon;
        let coversize_h = coversize as HighsInt;

        /* define the lifting function */
        let mut g = |z: f64| -> f64 {
            let hfrac = z / abar;
            let mut coef = 0.0f64;

            let mut h = (hfrac + 0.5).floor() as HighsInt;
            if h != 0 && (hfrac - h as f64).abs() * cpp_max(1.0, abar) <= epsilon && h < cplussize {
                halfintegral = true;
                coef = 0.5;
            }

            h = cpp_max(h - 1, 0);
            while h < coversize_h {
                if z <= s[h as usize] + feastol {
                    break;
                }
                h += 1;
            }

            coef + h as f64
        };

        self.rhs = HighsCDouble::from((coversize_h - 1) as f64);

        for i in 0..self.rowlen as usize {
            if vals[i] == 0.0 {
                continue;
            }
            if coverflag[i] == -1 {
                vals[i] = 1.0;
            } else {
                vals[i] = g(vals[i]);
            }
        }

        if halfintegral {
            self.rhs *= 2.0;
            for i in 0..self.rowlen as usize {
                vals[i] *= 2.0;
            }
        }

        // resulting cut is always integral
        self.integral_support = true;
        self.integral_coefficients = true;
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::separateLiftedMixedBinaryCover (lines 216-278)
    /// C++: `bool separateLiftedMixedBinaryCover()`
    #[allow(clippy::needless_range_loop)]
    pub fn separate_lifted_mixed_binary_cover(&mut self, vals: &mut [f64]) -> bool {
        // initialize indicators
        self.integral_support = false;
        self.integral_coefficients = false;

        let coversize = self.cover.len();
        let mut s = vec![0.0f64; coversize];
        let mut coverflag = vec![0u8; self.rowlen as usize];

        if coversize == 0 {
            return false;
        }

        for i in 0..coversize {
            coverflag[self.cover[i] as usize] = 1;
        }

        pdqsort_branchless_by(&mut self.cover[..], |&a, &b| {
            vals[a as usize] > vals[b as usize]
        });
        let mut sum = HighsCDouble::from(0.0);

        let lambda = self.lambda;
        let epsilon = self.epsilon;
        let mut p = coversize;
        for i in 0..coversize {
            let vc = vals[self.cover[i] as usize];
            if vc - lambda <= epsilon {
                p = i;
                break;
            }
            sum += vc;
            s[i] = f64::from(sum);
        }
        if p == 0 {
            return false;
        }

        /* define the lifting function */
        let phi = |a: f64| -> f64 {
            for i in 0..p {
                if a <= s[i] - lambda {
                    return f64::from((i as f64) * lambda);
                }

                if a <= s[i] {
                    return f64::from(((i + 1) as f64) * lambda + (HighsCDouble::from(a) - s[i]));
                }
            }

            f64::from((p as f64) * lambda + (HighsCDouble::from(a) - s[p - 1]))
        };

        self.rhs = -lambda;

        self.integral_coefficients = false;
        self.integral_support = true;
        for i in 0..self.rowlen as usize {
            if self.isintegral[i] == 0 {
                if vals[i] < 0.0 {
                    self.integral_support = false;
                } else {
                    vals[i] = 0.0;
                }
                continue;
            }

            if coverflag[i] != 0 {
                vals[i] = cpp_min(vals[i], f64::from(lambda));
                self.rhs += vals[i];
            } else {
                vals[i] = phi(vals[i]);
            }
        }

        true
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::separateLiftedMixedIntegerCover (lines 280-500)
    /// C++: `bool separateLiftedMixedIntegerCover()`
    #[allow(clippy::needless_range_loop)]
    pub fn separate_lifted_mixed_integer_cover(&mut self, vals: &mut [f64]) -> bool {
        // initialize indicators
        self.integral_support = false;
        self.integral_coefficients = false;

        let coversize = self.cover.len();

        let mut coverflag = vec![0u8; self.rowlen as usize];
        for &i in self.cover.iter() {
            coverflag[i as usize] = 1;
        }

        pdqsort_branchless_by(&mut self.cover[..], |&a, &b| {
            vals[a as usize] > vals[b as usize]
        });

        let mut a: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); coversize];
        let mut u: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); coversize + 1];
        let mut m: Vec<HighsCDouble> = vec![HighsCDouble::from(0.0); coversize + 1];

        let mut usum = HighsCDouble::from(0.0);
        let mut msum = HighsCDouble::from(0.0);
        // set up the partial sums of the upper bounds, and the contributions
        for cc in 0..coversize {
            let i = self.cover[cc] as usize;

            u[cc] = usum;
            m[cc] = msum;
            a[cc] = HighsCDouble::from(vals[i]);
            let ub = self.upper[i];
            usum += ub;
            msum += ub * a[cc];
        }

        u[coversize] = usum;
        m[coversize] = msum;

        let feastol = self.feastol;
        let lambda = self.lambda;

        // determine which variable in the cover we want to create the MIR inequality
        // from which we lift we try to select a variable to have the highest chance
        // of satisfying the facet conditions for the superadditive lifting function
        // gamma to be satisfied.
        let mut lpos: Option<usize> = None;
        let mut best_l_cplusend: usize = 0;
        let mut best_l_val: f64 = 0.0;
        let mut best_l_at_upper = true;

        for i in 0..coversize {
            let j = self.cover[i] as usize;
            let ub = self.upper[j];

            let at_upper = self.solval[j] >= ub - feastol;
            if at_upper && !best_l_at_upper {
                continue;
            }

            let mju = ub * vals[j];
            let mu = mju - lambda;

            if mu <= 10.0 * feastol {
                continue;
            }
            if vals[j].abs() < 1000.0 * feastol {
                continue;
            }

            let mudival = f64::from(mu / vals[j]);
            if HighsIntegers::is_integral(mudival, feastol) {
                continue;
            }
            let eta = mudival.ceil();

            let ulminusetaplusone = HighsCDouble::from(ub) - eta + 1.0;
            let cplusthreshold = ulminusetaplusone * vals[j];

            // std::upper_bound over the (sorted) cover with `value > vals[i]`
            let cthr = f64::from(cplusthreshold);
            let cplusend = self.cover.partition_point(|&e| !(cthr > vals[e as usize]));

            let mut mcplus = m[cplusend];
            if i < cplusend {
                mcplus -= mju;
            }

            let jl_val = f64::from(mcplus + eta * vals[j]);

            if jl_val > best_l_val || (!at_upper && best_l_at_upper) {
                lpos = Some(i);
                best_l_cplusend = cplusend;
                best_l_val = jl_val;
                best_l_at_upper = at_upper;
            }
        }

        let lpos = match lpos {
            Some(p) => p,
            None => return false,
        };

        let l = self.cover[lpos] as usize;
        let al = HighsCDouble::from(vals[l]);
        let upperl = self.upper[l];
        let mlu = upperl * al;
        let mu = mlu - lambda;

        a.truncate(best_l_cplusend);
        self.cover.truncate(best_l_cplusend);
        u.truncate(best_l_cplusend + 1);
        m.truncate(best_l_cplusend + 1);

        if lpos < best_l_cplusend {
            a.remove(lpos);
            self.cover.remove(lpos);
            u.remove(lpos + 1);
            m.remove(lpos + 1);
            for i in (lpos + 1)..best_l_cplusend {
                u[i] -= upperl;
                m[i] -= mlu;
            }
        }

        let cplussize = a.len();

        let mudival = f64::from(mu / al);
        let eta = mudival.ceil();
        let mut r = mu - mudival.floor() * al;
        // we multiply with r and it is important that it does not flip the sign
        // so we safe guard against tiny numerical errors here
        if r < 0.0 {
            r = HighsCDouble::from(0.0);
        }

        let ulminusetaplusone = HighsCDouble::from(upperl) - eta + 1.0;
        let cplusthreshold = ulminusetaplusone * al;

        let kmin = (eta - upperl - 0.5).floor() as HighsInt;
        let al_f = f64::from(al);

        let phi_l = |a_v: f64| -> f64 {
            let mut k: i64 = ((a_v / al_f) as i64).min(-1);

            while k >= kmin as i64 {
                if a_v >= (k as f64) * al + r {
                    return f64::from(a_v - ((k + 1) as f64) * r);
                }

                if a_v >= (k as f64) * al {
                    return f64::from((k as f64) * (al - r));
                }
                k -= 1;
            }

            f64::from((kmin as f64) * (al - r))
        };

        let kmax = (upperl - eta + 0.5).floor() as i64;

        let upper = &self.upper;
        let cover = &self.cover;
        let gamma_l = |z: f64| -> f64 {
            for i in 0..cplussize {
                let upperi = upper[cover[i] as usize] as HighsInt;

                for h in 0..=upperi {
                    let mih = m[i] + (h as f64) * a[i];
                    let uih = u[i] + (h as f64);
                    let mihplusdeltai = mih + a[i] - cplusthreshold;
                    if z <= mihplusdeltai {
                        return f64::from(uih * ulminusetaplusone * (al - r));
                    }

                    let mut k: i64 = f64::from((z - mihplusdeltai) / al) as i64 - 1;
                    while k <= kmax {
                        if z <= mihplusdeltai + (k as f64) * al + r {
                            return f64::from((uih * ulminusetaplusone + (k as f64)) * (al - r));
                        }

                        if z <= mihplusdeltai + ((k + 1) as f64) * al {
                            return f64::from(
                                (uih * ulminusetaplusone) * (al - r) + z - mih - a[i]
                                    + cplusthreshold
                                    - ((k + 1) as f64) * r,
                            );
                        }
                        k += 1;
                    }
                }
            }

            let mut p: i64 = f64::from((z - m[cplussize]) / al) as i64 - 1;
            loop {
                if z <= m[cplussize] + (p as f64) * al + r {
                    return f64::from((u[cplussize] * ulminusetaplusone + (p as f64)) * (al - r));
                }

                if z <= m[cplussize] + ((p + 1) as f64) * al {
                    return f64::from(
                        (u[cplussize] * ulminusetaplusone) * (al - r) + z
                            - m[cplussize]
                            - ((p + 1) as f64) * r,
                    );
                }
                p += 1;
            }
        };

        self.rhs = (HighsCDouble::from(upperl) - eta) * r - lambda;
        self.integral_support = true;
        self.integral_coefficients = false;
        for i in 0..self.rowlen as usize {
            if vals[i] == 0.0 {
                continue;
            }
            if self.isintegral[i] == 0 {
                if vals[i] < 0.0 {
                    self.integral_support = false;
                } else {
                    vals[i] = 0.0;
                }
                continue;
            }

            if coverflag[i] != 0 {
                vals[i] = -phi_l(-vals[i]);
                self.rhs += vals[i] * self.upper[i];
            } else {
                vals[i] = gamma_l(vals[i]);
            }
        }

        true
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::cmirCutGenerationHeuristic (lines 504-739)
    /// C++: `bool cmirCutGenerationHeuristic(double minEfficacy, bool onlyInitialCMIRScale = false)`
    /// C++ default arguments: `only_initial_cmir_scale = false`
    #[allow(clippy::needless_range_loop)]
    pub fn cmir_cut_generation_heuristic(
        &mut self,
        vals: &mut [f64],
        min_efficacy: f64,
        only_initial_cmir_scale: bool,
    ) -> bool {
        const MAX_CMIR_SCALE: f64 = 1e6;
        const F0MIN: f64 = 0.005;
        const F0MAX: f64 = 0.995;

        // initialize indicators
        self.integral_support = false;
        self.integral_coefficients = false;

        let mut continuouscontribution = 0.0f64;
        let mut continuoussqrnorm = 0.0f64;

        self.deltas.clear();
        self.integerinds.clear();

        let mut maxabsdelta = 0.0f64;

        self.complementation.resize(self.rowlen as usize, 0);

        for i in 0..self.rowlen {
            let iu = i as usize;
            if self.isintegral[iu] != 0 {
                self.integerinds.push(i);

                if self.upper[iu] < 2.0 * self.solval[iu] {
                    self.flip_complementation(vals, i);
                }

                if only_initial_cmir_scale {
                    continue;
                }

                if self.solval[iu] > self.feastol {
                    let delta = vals[iu].abs();
                    if delta <= 1e-4 || delta == maxabsdelta {
                        continue;
                    }
                    maxabsdelta = cpp_max(maxabsdelta, delta);
                    self.deltas.push(delta);
                }
            } else {
                let aj = vals[iu];
                self.update_violation_and_norm(
                    i,
                    aj,
                    &mut continuouscontribution,
                    &mut continuoussqrnorm,
                );
            }
        }

        if continuoussqrnorm == 0.0 && self.deltas.len() > 1 {
            let int_scale =
                HighsIntegers::integral_scale_vec(&self.deltas, self.feastol, kHighsTiny);

            if int_scale != 0.0 && int_scale <= 1e4 {
                let scalrhs = f64::from(self.rhs) * int_scale;
                let downrhs = fast_floor(scalrhs);

                let f0 = scalrhs - downrhs;
                if f0 >= F0MIN && f0 <= F0MAX {
                    self.deltas.push(1.0 / int_scale);
                }
            }
        }

        self.deltas.push(cpp_min(1.0, self.initial_scale));
        if !only_initial_cmir_scale {
            self.deltas
                .push(maxabsdelta + cpp_min(1.0, self.initial_scale));
        }

        // gather the per-integer data read by the scale loops below (vals and
        // solval/upper are not modified until the flip loop)
        self.cmir_gather.clear();
        for &j in self.integerinds.iter() {
            let ju = j as usize;
            self.cmir_gather
                .push((vals[ju], self.solval[ju], self.upper[ju]));
        }

        if self.deltas.is_empty() {
            return false;
        }

        pdqsort_by(&mut self.deltas[..], |a: &f64, b: &f64| a < b);
        let mut curdelta = self.deltas[0];
        for i in 1..self.deltas.len() {
            if self.deltas[i] - curdelta <= 10.0 * self.feastol {
                self.deltas[i] = 0.0;
            } else {
                curdelta = self.deltas[i];
            }
        }

        self.deltas.retain(|&d| d != 0.0);

        let mut bestdelta: f64 = -1.0;
        let mut bestefficacy = min_efficacy;

        for &delta in self.deltas.iter() {
            let scale = 1.0 / delta;
            let scalrhs = f64::from(self.rhs) * scale;
            let downrhs = fast_floor(scalrhs);

            let f0 = scalrhs - downrhs;
            if f0 < F0MIN || f0 > F0MAX {
                continue;
            }
            let oneoveroneminusf0 = 1.0 / (1.0 - f0);
            if oneoveroneminusf0 > MAX_CMIR_SCALE {
                continue;
            }

            let contscale = scale * oneoveroneminusf0;
            let mut sqrnorm = contscale * contscale * continuoussqrnorm;
            let mut viol = contscale * continuouscontribution - downrhs;

            for &(valj, solvalj, upperj) in self.cmir_gather.iter() {
                let scalaj = valj * scale;
                let downaj = fast_floor(scalaj + kHighsTiny);
                let fj = scalaj - downaj;
                let aj = downaj + cpp_max(0.0, (fj - f0) * oneoveroneminusf0);
                violation_and_norm_term(aj, solvalj, upperj, self.feastol, &mut viol, &mut sqrnorm);
            }

            // TODO: Add some minimum violation check, e.g., 0.001 * feastol?
            let efficacy = viol / sqrnorm.sqrt();
            if efficacy > bestefficacy {
                bestdelta = delta;
                bestefficacy = efficacy;
            }
        }

        if bestdelta == -1.0 {
            return false;
        }

        /* try if multiplying best delta by 2 4 or 8 gives a better efficacy */
        for k in 1..=3 {
            if only_initial_cmir_scale {
                break;
            }
            let delta = bestdelta * ((1 << k) as f64);
            let scale = 1.0 / delta;
            let scalrhs = f64::from(self.rhs) * scale;
            let downrhs = fast_floor(scalrhs);
            let f0 = scalrhs - downrhs;
            if f0 < F0MIN || f0 > F0MAX {
                continue;
            }

            let oneoveroneminusf0 = 1.0 / (1.0 - f0);
            if oneoveroneminusf0 > MAX_CMIR_SCALE {
                continue;
            }

            let contscale = scale * oneoveroneminusf0;
            let mut sqrnorm = contscale * contscale * continuoussqrnorm;
            let mut viol = contscale * continuouscontribution - downrhs;

            for &(valj, solvalj, upperj) in self.cmir_gather.iter() {
                let scalaj = valj * scale;
                let downaj = fast_floor(scalaj + kHighsTiny);
                let fj = scalaj - downaj;
                let aj = downaj + cpp_max(0.0, (fj - f0) * oneoveroneminusf0);
                violation_and_norm_term(aj, solvalj, upperj, self.feastol, &mut viol, &mut sqrnorm);
            }

            let efficacy = viol / sqrnorm.sqrt();
            if efficacy > bestefficacy {
                bestdelta = delta;
                bestefficacy = efficacy;
            }
        }

        // try to flip complementation of integers to increase efficacy
        for idx in 0..self.integerinds.len() {
            let k = self.integerinds[idx];
            let ku = k as usize;
            if self.upper[ku] == kHighsInf {
                continue;
            }
            if self.solval[ku] <= self.feastol {
                continue;
            }

            self.flip_complementation(vals, k);

            let delta = bestdelta;
            let scale = 1.0 / delta;
            let scalrhs = f64::from(self.rhs) * scale;
            let downrhs = fast_floor(scalrhs);

            let f0 = scalrhs - downrhs;
            if f0 < F0MIN || f0 > F0MAX {
                self.flip_complementation(vals, k);
                continue;
            }

            let oneoveroneminusf0 = 1.0 / (1.0 - f0);
            if oneoveroneminusf0 > MAX_CMIR_SCALE {
                self.flip_complementation(vals, k);
                continue;
            }

            let contscale = scale * oneoveroneminusf0;
            let mut sqrnorm = contscale * contscale * continuoussqrnorm;
            let mut viol = contscale * continuouscontribution - downrhs;

            for jdx in 0..self.integerinds.len() {
                let j = self.integerinds[jdx];
                let scalaj = vals[j as usize] * scale;
                let downaj = fast_floor(scalaj + kHighsTiny);
                let fj = scalaj - downaj;
                let aj = downaj + cpp_max(0.0, (fj - f0) * oneoveroneminusf0);
                self.update_violation_and_norm(j, aj, &mut viol, &mut sqrnorm);
            }

            let efficacy = viol / sqrnorm.sqrt();
            if efficacy > bestefficacy {
                bestefficacy = efficacy;
            } else {
                self.flip_complementation(vals, k);
            }
        }

        let scale = 1.0 / HighsCDouble::from(bestdelta);
        let scalrhs = self.rhs * scale;
        let downrhs = f64::from(scalrhs).floor();

        let f0 = scalrhs - downrhs;
        let oneoveroneminusf0 = 1.0 / (1.0 - f0);

        self.rhs = HighsCDouble::from(downrhs * bestdelta);
        self.integral_support = true;
        self.integral_coefficients = false;
        for j in 0..self.rowlen as usize {
            if vals[j] == 0.0 {
                continue;
            }
            if self.isintegral[j] == 0 {
                if vals[j] > 0.0 {
                    vals[j] = 0.0;
                } else {
                    vals[j] = f64::from(vals[j] * oneoveroneminusf0);
                    self.integral_support = false;
                }
            } else {
                let scalaj = scale * vals[j];
                let downaj = f64::from(scalaj + kHighsTiny).floor();
                let fj = scalaj - downaj;
                let mut aj = HighsCDouble::from(downaj);
                if fj > f0 {
                    aj += (fj - f0) * oneoveroneminusf0;
                }

                vals[j] = f64::from(aj * bestdelta);
            }
        }

        // HiGHS: `#ifndef NDEBUG` check that the computed cut has the correct efficacy
        if cfg!(debug_assertions) {
            let mut checkviol = -downrhs * bestdelta;
            let mut checknorm = 0.0f64;
            for j in 0..self.rowlen as usize {
                if vals[j] == 0.0 {
                    continue;
                }
                let vj = vals[j];
                self.update_violation_and_norm(j as HighsInt, vj, &mut checkviol, &mut checknorm);
            }
            if checknorm != 0.0 {
                // the efficacy can become infinite if the cut 0 <= -1 is derived
                let checkefficacy = checkviol / checknorm.sqrt();
                debug_assert!(
                    (checkefficacy - bestefficacy).abs() < 0.001 || checkefficacy.abs() >= 1e30
                );
            }
        }

        true
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::scale (lines 741-758)
    /// C++: `double scale(double val)`
    pub fn scale(&mut self, vals: &mut [f64], val: f64) -> f64 {
        let (_, exp) = frexp(val);
        let mut expshift = -exp;

        // Don't scale the coefficients by more than +1024 (violations can increase)
        expshift = cpp_min(10, expshift);

        // Scale rhs
        self.rhs = self.rhs.ldexp(expshift);

        // Scale row
        for i in 0..self.rowlen as usize {
            vals[i] = ldexp_f64(vals[i], expshift);
        }

        // Return scaling factor
        ldexp_f64(1.0, expshift)
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::postprocessCut (lines 760-899)
    /// C++: `bool postprocessCut(const HighsDomain & globaldom)`
    // PORT NOTE: `globaldom` is the caller's domain (`ms.mipdata_.domains[..]`, passed
    // attached by the caller). `lpRelaxation.isColIntegral` takes `ms`.
    #[allow(clippy::needless_range_loop)]
    pub fn postprocess_cut(
        &mut self,
        ms: &HighsMipSolver,
        lp_relaxation: &HighsLpRelaxation,
        vals: &mut [f64],
        inds: &mut [HighsInt],
        globaldom: &HighsDomain,
    ) -> bool {
        // right hand sides slightly below zero are likely due to numerical errors and
        // can cause numerical troubles with scaling, so set them to zero
        if self.rhs < 0.0 && self.rhs > -self.epsilon {
            self.rhs = HighsCDouble::from(0.0);
        }

        if self.integral_support && self.integral_coefficients {
            // if the cut is known to be integral no postprocessing is needed and we
            // simply remove zero coefficients
            for i in (0..self.rowlen).rev() {
                let iu = i as usize;
                if vals[iu] == 0.0 {
                    self.rowlen -= 1;
                    let rl = self.rowlen as usize;
                    inds[iu] = inds[rl];
                    vals[iu] = vals[rl];
                }
            }
            return true;
        }

        // determine maximal absolute coefficient
        let mut max_abs_value = 0.0f64;
        for i in 0..self.rowlen as usize {
            max_abs_value = cpp_max(vals[i].abs(), max_abs_value);
        }

        // determine minimal allowed coefficient
        let min_coefficient_value = 100.0 * self.feastol * cpp_max(max_abs_value, 1e-3);

        // remove small coefficients and check whether the remaining support is
        // integral
        self.integral_support = true;
        for i in (0..self.rowlen).rev() {
            let iu = i as usize;
            if vals[iu] == 0.0 {
                continue;
            }
            if vals[iu].abs() <= min_coefficient_value {
                if vals[iu] < 0.0 {
                    let ub = globaldom.col_upper_[inds[iu] as usize];
                    if ub == kHighsInf {
                        return false;
                    } else {
                        self.rhs -= ub * vals[iu];
                    }
                } else {
                    let lb = globaldom.col_lower_[inds[iu] as usize];
                    if lb == -kHighsInf {
                        return false;
                    } else {
                        self.rhs -= lb * vals[iu];
                    }
                }

                vals[iu] = 0.0;
                continue;
            }

            if self.integral_support && !lp_relaxation.is_col_integral(ms, inds[iu]) {
                self.integral_support = false;
            }
        }

        // remove zeros in place
        for i in (0..self.rowlen).rev() {
            let iu = i as usize;
            if vals[iu] == 0.0 {
                self.rowlen -= 1;
                let rl = self.rowlen as usize;
                inds[iu] = inds[rl];
                vals[iu] = vals[rl];
            }
        }

        if self.rowlen == 0 {
            return false;
        }

        if self.integral_support {
            // integral support -> determine scale to make all coefficients integral
            let intscale =
                HighsIntegers::integral_scale(vals, self.rowlen, self.feastol, self.epsilon);

            let mut scale_smallest_val_to_one = true;

            if intscale != 0.0 && intscale * cpp_max(1.0, max_abs_value) <= (1u64 << 52) as f64 {
                // A scale to make all value integral was found. The scale is only
                // rejected if it is in a range where not all integral values are
                // representable in double precision anymore or cannot be correctly
                // rounded by adding 0.5 and casting to an int. The latter starts at 2^52
                // + 1 since adding 0.5 will round upwards to the next even number for
                // that magnitude. So the largest acceptable value is 2^52 and when at
                // most that value we want to always use the scale to adjust the
                // coefficients and right hand side for numerical safety reasons. If the
                // resulting integral values are too large, however, we scale the cut down
                // by shifting the exponent.
                self.rhs.renormalize();
                self.rhs *= intscale;
                max_abs_value = HighsIntegers::nearest_integer(max_abs_value * intscale) as f64;
                for i in 0..self.rowlen as usize {
                    let scaleval = intscale * HighsCDouble::from(vals[i]);
                    let intval = HighsIntegers::nearest_integer(f64::from(scaleval)) as f64;
                    let delta = f64::from(scaleval - intval);

                    vals[i] = intval;

                    // if the coefficient would be strengthened by rounding, we add the
                    // upperbound constraint to make it exactly integral instead and
                    // therefore weaken the right hand side
                    if delta < 0.0 {
                        let ub = globaldom.col_upper_[inds[i] as usize];
                        if ub == kHighsInf {
                            return false;
                        }

                        self.rhs -= delta * ub;
                    } else {
                        let lb = globaldom.col_lower_[inds[i] as usize];
                        if lb == -kHighsInf {
                            return false;
                        }

                        self.rhs -= delta * lb;
                    }
                }

                // finally we can round down the right hand side. Therefore in most cases
                // small errors for which the upper bound constraints where used and the
                // right hand side was weakened, do not weaken the final cut.
                self.rhs = (self.rhs + self.feastol).floor();

                if intscale * max_abs_value * self.feastol < 0.5 {
                    // integral scale leads to small enough values, accept scale
                    scale_smallest_val_to_one = false;
                    self.integral_coefficients = true;
                }
            }

            if scale_smallest_val_to_one {
                // integral scale lead to very large coefficient values. We now shift the
                // exactly integral values down such that the smallest coefficient is
                // around 1
                let mut min_abs_value = kHighsInf;
                for i in 0..self.rowlen as usize {
                    min_abs_value = cpp_min(vals[i].abs(), min_abs_value);
                }

                self.scale(vals, min_abs_value - self.epsilon);
            }
        } else {
            // the support is not integral, scale cut to have the largest coefficient
            // around 1.0
            self.scale(vals, max_abs_value - self.epsilon);
        }

        true
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::preprocessBaseInequality (lines 901-1037)
    /// HiGHS 1.16 (#3267): HighsCutGeneration.cpp::HighsCutGeneration::preprocessBaseInequality (lines 901-1040)
    /// C++: `bool preprocessBaseInequality(bool & hasUnboundedInts, bool & hasGeneralInts, bool & hasContinuous)`
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::needless_range_loop)]
    pub fn preprocess_base_inequality(
        &mut self,
        ms: &HighsMipSolver,
        lp_relaxation: &HighsLpRelaxation,
        vals: &mut [f64],
        inds: &mut [HighsInt],
        has_unbounded_ints: &mut bool,
        has_general_ints: &mut bool,
        has_continuous: &mut bool,
    ) -> bool {
        // preprocess the inequality before cut generation
        // 1. Determine the maximal activity to check for trivial redundancy and
        // tighten coefficients
        // 2. Check for presence of continuous variables and unbounded integers as not
        // all methods for cut generation are applicable in that case
        // 3. Remove coefficients that are below the feasibility tolerance to avoid
        // numerical troubles, use bound constraints to cancel them and
        // reject base inequalities where that is not possible due to unbounded
        // variables
        *has_unbounded_ints = false;
        *has_continuous = false;
        *has_general_ints = false;
        let mut num_zeros: HighsInt = 0;

        let mut maxact = -self.feastol;
        let mut max_abs_val = 0.0f64;
        for i in 0..self.rowlen as usize {
            max_abs_val = cpp_max(vals[i].abs(), max_abs_val);
        }

        self.initial_scale = self.scale(vals, max_abs_val);

        let rowlen = self.rowlen;
        self.isintegral.resize(rowlen as usize, 0);
        for i in 0..rowlen {
            let iu = i as usize;
            // we do not want to have integral variables with small coefficients as this
            // may lead to numerical instabilities during cut generation
            // Therefore we relax integral variables with small coefficients to
            // continuous ones because they still might have a non-negligible
            // contribution e.g. when they come from integral rows with coefficients on
            // the larger side. When we relax them to continuous variables they will be
            // complemented so that their solution value is closest to zero and then
            // will be relaxed if their value is positive or their maximal contribution
            // is below feasibility tolerance.
            self.isintegral[iu] = (lp_relaxation.is_col_integral(ms, inds[iu])
                && vals[iu].abs() > 10.0 * self.feastol) as u8;

            if self.isintegral[iu] == 0 {
                // complement non-integer variable (cmir separation heuristic complements
                // integral variables in the same way)
                if self.upper[iu] < 2.0 * self.solval[iu] {
                    if self.complementation.is_empty() {
                        self.complementation.resize(rowlen as usize, 0);
                    }

                    self.flip_complementation(vals, i);
                }

                // relax positive continuous variables and those with small contributions
                if vals[iu] > 0.0 || vals[iu].abs() * self.upper[iu] <= 10.0 * self.feastol {
                    if vals[iu] < 0.0 {
                        if self.upper[iu] == kHighsInf {
                            return false;
                        }
                        self.rhs -= vals[iu] * self.upper[iu];
                    }

                    num_zeros += 1;
                    vals[iu] = 0.0;
                    continue;
                }

                *has_continuous = true;
            } else {
                // HiGHS 1.16 #3267 (h16_big_act): HighsCutGeneration.cpp lines 963-965 and 966 (1.16.0)
                let h16_big_act = ms.options_mip_.h16_big_act;
                if h16_big_act
                    && self.upper[iu] != kHighsInf
                    && (vals[iu] * self.upper[iu]).abs() > kHighsInvMacheps
                {
                    return false;
                }
                let unbounded_int = if h16_big_act {
                    self.upper[iu] >= kHighsIInf as f64
                } else {
                    self.upper[iu] == kHighsInf
                };
                if unbounded_int {
                    *has_unbounded_ints = true;
                    *has_general_ints = true;
                } else if self.upper[iu] != 1.0 {
                    *has_general_ints = true;
                }

                if vals[iu] > 0.0 {
                    maxact += vals[iu] * self.upper[iu];
                }
            }
        }

        let max_len: HighsInt = if self.loosen_len {
            HighsInt::MAX
        } else {
            (100.0 + 0.15 * (lp_relaxation.num_cols() as f64)) as HighsInt
        };

        if rowlen - num_zeros > max_len {
            let num_cancel: HighsInt = rowlen - num_zeros - max_len;
            let mut cancel_nzs: Vec<HighsInt> = Vec::new();

            for i in 0..rowlen as usize {
                let cancel_slack = if vals[i] > 0.0 {
                    self.solval[i]
                } else {
                    self.upper[i] - self.solval[i]
                };
                if cancel_slack <= self.feastol {
                    cancel_nzs.push(i as HighsInt);
                }
            }

            if (cancel_nzs.len() as HighsInt) < num_cancel {
                return false;
            }
            if (cancel_nzs.len() as HighsInt) > num_cancel {
                cpp_partial_sort_by(&mut cancel_nzs[..], num_cancel as usize, |&a, &b| {
                    vals[a as usize].abs() < vals[b as usize].abs()
                });
            }

            for i in 0..num_cancel as usize {
                let j = cancel_nzs[i] as usize;

                if vals[j] < 0.0 {
                    self.rhs -= vals[j] * self.upper[j];
                } else {
                    maxact -= vals[j] * self.upper[j];
                }

                vals[j] = 0.0;
            }

            num_zeros += num_cancel;
        }

        if num_zeros != 0 {
            // remove zeros in place
            // (C++ has two loops, one without and one with the complementation copy; the
            // complementation vector cannot become empty inside the loop, so one loop with
            // a guarded copy is the same)
            let has_complementation = !self.complementation.is_empty();
            for i in (0..self.rowlen).rev() {
                let iu = i as usize;
                if vals[iu] == 0.0 {
                    self.rowlen -= 1;
                    let rl = self.rowlen as usize;
                    inds[iu] = inds[rl];
                    vals[iu] = vals[rl];
                    self.upper[iu] = self.upper[rl];
                    self.solval[iu] = self.solval[rl];
                    self.isintegral[iu] = self.isintegral[rl];
                    if has_complementation {
                        self.complementation[iu] = self.complementation[rl];
                    }
                    num_zeros -= 1;
                    if num_zeros == 0 {
                        break;
                    }
                }
            }
        }

        maxact > self.rhs
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::generateCut (lines 1058-1197)
    /// C++: `bool generateCut(HighsTransformedLp & transLp, std::vector<HighsInt> & inds_, std::vector<double> & vals_, double & rhs_, bool onlyInitialCMIRScale = false)`
    /// C++ default arguments: `only_initial_cmir_scale = false`
    // PORT NOTE: `ms` is `&mut` for `transLp.transform` / `untransform` (they detach the
    // implications) and `addCut`. `transLp.getGlobaldom().tightenCoefficients(..)` is
    // `trans_lp.get_globaldom(ms).tighten_coefficients(ms, ..)`.
    #[allow(clippy::too_many_arguments)]
    pub fn generate_cut(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &HighsLpRelaxation,
        trans_lp: &mut HighsTransformedLp,
        inds_: &mut Vec<HighsInt>,
        vals_: &mut Vec<f64>,
        rhs_: &mut f64,
        only_initial_cmir_scale: bool,
    ) -> bool {
        let mut ints_positive = true;
        if !trans_lp.transform(
            ms,
            lp_relaxation,
            vals_,
            &mut self.upper,
            &mut self.solval,
            inds_,
            rhs_,
            &mut ints_positive,
            false,
            false,
        ) {
            return false;
        }

        self.rowlen = inds_.len() as HighsInt;
        self.rhs = HighsCDouble::from(*rhs_);
        self.complementation.clear();
        let mut has_unbounded_ints = false;
        let mut has_general_ints = false;
        let mut has_continuous = false;
        if !self.preprocess_base_inequality(
            ms,
            lp_relaxation,
            &mut vals_[..],
            &mut inds_[..],
            &mut has_unbounded_ints,
            &mut has_general_ints,
            &mut has_continuous,
        ) {
            return false;
        }

        // it can happen that there is an unbounded integer variable during the
        // transform call so that the integers are not transformed to positive values.
        // Now the call to preprocessBaseInequality may have removed the unbounded
        // integer, e.g. due to a small coefficient value, so that we can still use
        // the lifted inequalities instead of cmir. We need to make sure, however,
        // that the cut values are transformed to positive coefficients first, which
        // we do below.
        if !has_unbounded_ints && !ints_positive {
            self.complementation.resize(self.rowlen as usize, 0);

            for i in 0..self.rowlen as usize {
                if vals_[i] > 0.0 || self.isintegral[i] == 0 {
                    continue;
                }

                self.flip_complementation(&mut vals_[..], i as HighsInt);
            }
        }

        // try to generate a cut
        if !self.try_generate_cut(
            ms,
            inds_,
            vals_,
            has_unbounded_ints,
            has_general_ints,
            has_continuous,
            10.0 * self.feastol,
            only_initial_cmir_scale,
            true,
            true,
        ) {
            return false;
        }

        // remove the complementation if exists
        self.remove_complementation(&mut vals_[..]);

        // remove zeros in place
        for i in (0..self.rowlen).rev() {
            let iu = i as usize;
            if vals_[iu] == 0.0 {
                self.rowlen -= 1;
                let rl = self.rowlen as usize;
                inds_[iu] = inds_[rl];
                vals_[iu] = vals_[rl];
            }
        }

        // transform the cut back into the original space, i.e. remove the bound
        // substitution and replace implicit slack variables
        *rhs_ = f64::from(self.rhs);
        vals_.resize(self.rowlen as usize, 0.0);
        inds_.resize(self.rowlen as usize, 0);
        if !trans_lp.untransform(ms, lp_relaxation, vals_, inds_, rhs_, false) {
            return false;
        }

        self.rowlen = inds_.len() as HighsInt;
        self.rhs = HighsCDouble::from(*rhs_);

        ms.mipdata_
            .debug_solution
            .check_cut(&inds_[..], &vals_[..], self.rowlen, *rhs_);
        // apply cut postprocessing including scaling and removal of small
        // coefficients
        if !self.postprocess_cut(
            &*ms,
            lp_relaxation,
            &mut vals_[..],
            &mut inds_[..],
            trans_lp.get_globaldom(&*ms),
        ) {
            return false;
        }
        *rhs_ = f64::from(self.rhs);
        vals_.resize(self.rowlen as usize, 0.0);
        inds_.resize(self.rowlen as usize, 0);

        ms.mipdata_
            .debug_solution
            .check_cut(&inds_[..], &vals_[..], self.rowlen, *rhs_);

        // finally determine the violation of the cut in the original space
        let mut violation = HighsCDouble::from(-*rhs_);
        let sol = &lp_relaxation.get_solution().col_value;
        for i in 0..self.rowlen as usize {
            violation += sol[inds_[i] as usize] * vals_[i];
        }

        if violation <= 10.0 * self.feastol {
            return false;
        }

        trans_lp.get_globaldom(&*ms).tighten_coefficients(
            &*ms,
            &mut inds_[..],
            &mut vals_[..],
            self.rowlen,
            rhs_,
        );

        // if the cut is violated by a small factor above the feasibility
        // tolerance, add it to the cutpool
        let len = inds_.len() as HighsInt;
        let integral = self.integral_support && self.integral_coefficients;
        let cutindex = detach!(ms.mipdata_.cutpools[self.cutpool], |p| p.add_cut(
            ms,
            &mut inds_[..],
            &mut vals_[..],
            len,
            *rhs_,
            integral,
            true,
            true,
            false
        ));

        // only return true if cut was accepted by the cutpool, i.e. not a duplicate
        // of a cut already in the pool
        cutindex != -1
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::generateConflict (lines 1199-1291)
    /// C++: `bool generateConflict(const HighsDomain & localdomain, const HighsDomain & globaldom, std::vector<HighsInt> & proofinds, std::vector<double> & proofvals, double & proofrhs)`
    // PORT NOTE: `globaldom` is an index into `ms.mipdata_.domains` (read attached; the C++
    // passes `mipdata_->domain`); `localdomain` is the caller's detached node domain.
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::needless_range_loop)]
    pub fn generate_conflict(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &HighsLpRelaxation,
        localdomain: &HighsDomain,
        globaldom: usize,
        proofinds: &mut Vec<HighsInt>,
        proofvals: &mut Vec<f64>,
        proofrhs: &mut f64,
    ) -> bool {
        // rusolve ablation lever mip_conflict_off_mask bit 32
        if ms.options_mip_.mip_conflict_off_mask & 32 != 0 {
            return false;
        }
        self.rhs = HighsCDouble::from(*proofrhs);
        self.rowlen = proofinds.len() as HighsInt;

        ms.mipdata_.debug_solution.check_cut(
            &proofinds[..],
            &proofvals[..],
            self.rowlen,
            *proofrhs,
        );

        let rowlen = self.rowlen as usize;
        self.complementation = vec![0u8; rowlen];

        self.upper.resize(rowlen, 0.0);
        self.solval.resize(rowlen, 0.0);

        let gd: &HighsDomain = &ms.mipdata_.domains[globaldom];

        let mut activity = 0.0f64;
        for i in 0..rowlen {
            let col = proofinds[i] as usize;

            self.upper[i] = gd.col_upper_[col] - gd.col_lower_[col];

            self.solval[i] = if proofvals[i] < 0.0 {
                cpp_min(gd.col_upper_[col], localdomain.col_upper_[col])
            } else {
                cpp_max(gd.col_lower_[col], localdomain.col_lower_[col])
            };
            if proofvals[i] < 0.0 && gd.col_upper_[col] != kHighsInf {
                self.rhs -= gd.col_upper_[col] * proofvals[i];
                proofvals[i] = -proofvals[i];
                self.complementation[i] = 1;

                self.solval[i] = gd.col_upper_[col] - self.solval[i];
            } else {
                self.rhs -= gd.col_lower_[col] * proofvals[i];
                self.complementation[i] = 0;
                self.solval[i] -= gd.col_lower_[col];
            }

            activity += self.solval[i] * proofvals[i];
        }

        if activity > self.rhs {
            let sol_scale = f64::from(self.rhs) / activity;
            for i in 0..rowlen {
                self.solval[i] *= sol_scale;
            }
        }

        let mut has_unbounded_ints = false;
        let mut has_general_ints = false;
        let mut has_continuous = false;

        if !self.preprocess_base_inequality(
            &*ms,
            lp_relaxation,
            &mut proofvals[..],
            &mut proofinds[..],
            &mut has_unbounded_ints,
            &mut has_general_ints,
            &mut has_continuous,
        ) {
            return false;
        }

        // try to generate a cut
        if !self.try_generate_cut(
            &*ms,
            proofinds,
            proofvals,
            has_unbounded_ints,
            has_general_ints,
            has_continuous,
            self.feastol,
            false,
            false,
            false,
        ) {
            return false;
        }

        // remove the complementation
        if !self.complementation.is_empty() {
            for i in 0..self.rowlen as usize {
                if self.complementation[i] != 0 {
                    self.rhs -= gd.col_upper_[proofinds[i] as usize] * proofvals[i];
                    proofvals[i] = -proofvals[i];
                } else {
                    self.rhs += gd.col_lower_[proofinds[i] as usize] * proofvals[i];
                }
            }
        }

        // apply cut postprocessing including scaling and removal of small
        // coefficients
        if !self.postprocess_cut(
            &*ms,
            lp_relaxation,
            &mut proofvals[..],
            &mut proofinds[..],
            gd,
        ) {
            return false;
        }

        let rl = self.rowlen;
        proofvals.resize(rl as usize, 0.0);
        proofinds.resize(rl as usize, 0);
        *proofrhs = f64::from(self.rhs);

        let cutintegral = self.integral_support && self.integral_coefficients;

        gd.tighten_coefficients(&*ms, &mut proofinds[..], &mut proofvals[..], rl, proofrhs);

        let cutindex = detach!(ms.mipdata_.cutpools[self.cutpool], |p| p.add_cut(
            ms,
            &mut proofinds[..],
            &mut proofvals[..],
            rl,
            *proofrhs,
            cutintegral,
            true,
            true,
            true
        ));

        // only return true if cut was accepted by the cutpool, i.e. not a duplicate
        // of a cut already in the pool
        if cutindex != -1 {
            crate::work::add_event(crate::work::Event::ConfCutsAdded, 1);
        }
        cutindex != -1
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::finalizeAndAddCut (lines 1293-1349)
    /// C++: `bool finalizeAndAddCut(const HighsDomain & globaldom, std::vector<HighsInt> & inds_, std::vector<double> & vals_, double & rhs_)`
    // PORT NOTE: `globaldom` is an index into `ms.mipdata_.domains` (see `generate_conflict`).
    #[allow(clippy::needless_range_loop)]
    pub fn finalize_and_add_cut(
        &mut self,
        ms: &mut HighsMipSolver,
        lp_relaxation: &HighsLpRelaxation,
        globaldom: usize,
        inds_: &mut Vec<HighsInt>,
        vals_: &mut Vec<f64>,
        rhs_: &mut f64,
    ) -> bool {
        self.complementation.clear();
        self.rowlen = inds_.len() as HighsInt;
        self.rhs = HighsCDouble::from(*rhs_);

        self.integral_support = true;
        self.integral_coefficients = false;
        // remove zeros in place
        for i in (0..self.rowlen).rev() {
            let iu = i as usize;
            if vals_[iu] == 0.0 {
                self.rowlen -= 1;
                let rl = self.rowlen as usize;
                inds_[iu] = inds_[rl];
                vals_[iu] = vals_[rl];
            } else {
                self.integral_support &= lp_relaxation.is_col_integral(&*ms, inds_[iu]);
            }
        }

        vals_.resize(self.rowlen as usize, 0.0);
        inds_.resize(self.rowlen as usize, 0);

        ms.mipdata_
            .debug_solution
            .check_cut(&inds_[..], &vals_[..], self.rowlen, *rhs_);

        let gd: &HighsDomain = &ms.mipdata_.domains[globaldom];

        // apply cut postprocessing including scaling and removal of small
        // coefficients
        if !self.postprocess_cut(&*ms, lp_relaxation, &mut vals_[..], &mut inds_[..], gd) {
            return false;
        }
        *rhs_ = f64::from(self.rhs);
        vals_.resize(self.rowlen as usize, 0.0);
        inds_.resize(self.rowlen as usize, 0);

        ms.mipdata_
            .debug_solution
            .check_cut(&inds_[..], &vals_[..], self.rowlen, *rhs_);

        // finally determine the violation of the cut in the original space
        let mut violation = HighsCDouble::from(-*rhs_);
        let sol = &lp_relaxation.get_solution().col_value;
        for i in 0..self.rowlen as usize {
            violation += sol[inds_[i] as usize] * vals_[i];
        }

        if violation <= 10.0 * self.feastol {
            return false;
        }

        gd.tighten_coefficients(&*ms, &mut inds_[..], &mut vals_[..], self.rowlen, rhs_);

        // if the cut is violated by a small factor above the feasibility
        // tolerance, add it to the cutpool
        let len = inds_.len() as HighsInt;
        let integral = self.integral_support && self.integral_coefficients;
        let cutindex = detach!(ms.mipdata_.cutpools[self.cutpool], |p| p.add_cut(
            ms,
            &mut inds_[..],
            &mut vals_[..],
            len,
            *rhs_,
            integral,
            true,
            true,
            false
        ));

        // only return true if cut was accepted by the cutpool, i.e. not a duplicate
        // of a cut already in the pool
        cutindex != -1
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::flipComplementation (lines 1351-1360)
    /// C++: `void flipComplementation(HighsInt index)`
    pub fn flip_complementation(&mut self, vals: &mut [f64], index: HighsInt) {
        let i = index as usize;

        // only variables with finite upper bounds can be complemented
        // (C++ assert(upper[index] != kHighsInf) not mirrored)

        // flip complementation
        self.complementation[i] = 1 - self.complementation[i];
        self.solval[i] = self.upper[i] - self.solval[i];
        self.rhs -= self.upper[i] * vals[i];
        vals[i] = -vals[i];
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::removeComplementation (lines 1362-1368)
    /// C++: `void removeComplementation()`
    pub fn remove_complementation(&mut self, vals: &mut [f64]) {
        // remove the complementation if it exists
        if self.complementation.is_empty() {
            return;
        }
        for i in 0..self.rowlen {
            if self.complementation[i as usize] != 0 {
                self.flip_complementation(vals, i);
            }
        }
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::updateViolationAndNorm (lines 1370-1380)
    /// C++: `void updateViolationAndNorm(HighsInt index, double aj, double & violation, double & norm) const`
    pub fn update_violation_and_norm(
        &self,
        index: HighsInt,
        aj: f64,
        violation: &mut f64,
        norm: &mut f64,
    ) {
        let iu = index as usize;
        // update violation
        *violation += aj * self.solval[iu];

        // update norm
        if aj > 0.0 && self.solval[iu] <= self.feastol {
            return;
        }
        if aj < 0.0 && self.solval[iu] >= self.upper[iu] - self.feastol {
            return;
        }
        *norm += aj * aj;
    }

    /// HiGHS: HighsCutGeneration.cpp::HighsCutGeneration::tryGenerateCut (lines 1382-1482)
    /// C++: `bool tryGenerateCut(std::vector<HighsInt> & inds_, std::vector<double> & vals_, bool hasUnboundedInts, bool hasGeneralInts, bool hasContinuous, double minEfficacy, bool onlyInitialCMIRScale = false, bool allowRejectCut = true, bool lpSol = true)`
    /// C++ default arguments: `only_initial_cmir_scale = false`, `allow_reject_cut = true`, `lp_sol = true`
    // PORT NOTE: on entry the C++ `inds` / `vals` point at `inds_` / `vals_` (every caller set
    // them so). The copy into `tmpInds` / `tmpVals` and the re-point at them for the final
    // cmir call is: `mem::take` the two tmp vectors, run the cmir on `&mut tmp_vals[..]`, then
    // `swap` with `inds_` / `vals_` exactly where the C++ swaps, and put them back.
    #[allow(clippy::too_many_arguments)]
    pub fn try_generate_cut(
        &mut self,
        ms: &HighsMipSolver,
        inds_: &mut Vec<HighsInt>,
        vals_: &mut Vec<f64>,
        has_unbounded_ints: bool,
        has_general_ints: bool,
        has_continuous: bool,
        min_efficacy: f64,
        only_initial_cmir_scale: bool,
        allow_reject_cut: bool,
        lp_sol: bool,
    ) -> bool {
        // use cmir if there are unbounded integer variables
        if has_unbounded_ints {
            return self.cmir_cut_generation_heuristic(
                &mut vals_[..],
                min_efficacy,
                only_initial_cmir_scale,
            );
        }

        // 0. Save data before determining cover and applying lifting functions
        // C++ `tmpVals.assign(..)` / `tmpComplementation = complementation`: copy into the
        // existing buffers, so the capacity of the previous call is reused.
        let rl = self.rowlen as usize;
        self.tmp_vals.clear();
        self.tmp_vals.extend_from_slice(&vals_[..rl]);
        self.tmp_inds.clear();
        self.tmp_inds.extend_from_slice(&inds_[..rl]);
        self.tmp_complementation.clone_from(&self.complementation);
        self.tmp_solval.clone_from(&self.solval);
        let mut tmp_rhs = self.rhs;

        // 1. Determine a cover, cover does not need to be minimal as neither of
        //    the lifting functions have minimality of the cover as necessary facet
        //    condition
        let mut success = false;
        let mut save_integal_support = false;
        let mut save_integral_coefficients = false;
        // C++: do { ... } while (false), with `break` out of the block
        'cover: {
            if !self.determine_cover(ms, &mut vals_[..], &mut inds_[..], lp_sol) {
                break 'cover;
            }

            // 2. use superadditive lifting function depending on structure of base
            //    inequality:
            //    We have 3 lifting functions available for pure binary knapsack sets,
            //    for mixed-binary knapsack sets and for mixed integer knapsack sets.
            if !has_continuous && !has_general_ints {
                self.separate_lifted_knapsack_cover(&mut vals_[..]);
                success = true;
            } else if has_general_ints {
                success = self.separate_lifted_mixed_integer_cover(&mut vals_[..]);
            } else {
                success = self.separate_lifted_mixed_binary_cover(&mut vals_[..]);
            }
        }

        let mut min_mir_efficacy = min_efficacy;
        if success {
            // save data that might otherwise be overwritten when calling the cmir
            // separator
            save_integal_support = self.integral_support;
            save_integral_coefficients = self.integral_coefficients;

            // compute violation and squared norm
            let mut violation = -f64::from(self.rhs);
            let mut sqrnorm = 0.0f64;
            for i in 0..self.rowlen as usize {
                let aj = vals_[i];
                self.update_violation_and_norm(i as HighsInt, aj, &mut violation, &mut sqrnorm);
            }

            // compute efficacy (distance cut off)
            let efficacy = violation / sqrnorm.sqrt();
            if allow_reject_cut && efficacy <= min_efficacy {
                // reject cut
                success = false;
                self.rhs = tmp_rhs;
            } else {
                // accept cut and increase minimum efficiency requirement for cmir cut
                min_mir_efficacy += efficacy;
                std::mem::swap(&mut tmp_rhs, &mut self.rhs);
            }
        }

        // restore indices and values; lifting methods do not modify complementation
        // and, thus, complementation-related data does not have to be restored here.
        // The cmir runs on the saved (tmp) values.
        let mut tmp_vals = std::mem::take(&mut self.tmp_vals);

        if self.cmir_cut_generation_heuristic(
            &mut tmp_vals[..],
            min_mir_efficacy,
            only_initial_cmir_scale,
        ) {
            // take the cmir cut as it is better
            std::mem::swap(inds_, &mut self.tmp_inds);
            std::mem::swap(vals_, &mut tmp_vals);
            self.tmp_vals = tmp_vals;
            return true;
        }
        self.tmp_vals = tmp_vals;

        if success {
            // take the previous lifted cut as cmir could not improve
            // we restore the old complementation vector, right hand side, and values
            self.rhs = tmp_rhs;
            // note that the solution vector solval also needs to be restored because it
            // depends on the complementation. it would be OK not to restore solval, if
            // there would be a guarantee that it is not used from here on.
            std::mem::swap(&mut self.complementation, &mut self.tmp_complementation);
            std::mem::swap(&mut self.solval, &mut self.tmp_solval);
            // restore indicators
            self.integral_support = save_integal_support;
            self.integral_coefficients = save_integral_coefficients;
            return true;
        }

        // neither cmir nor lifted cut successful
        false
    }
}

/// The body of `HighsCutGeneration::updateViolationAndNorm` for a column whose
/// `solval` and `upper` are already in hand (same operations, same order).
#[inline(always)]
fn violation_and_norm_term(
    aj: f64,
    solval: f64,
    upper: f64,
    feastol: f64,
    violation: &mut f64,
    norm: &mut f64,
) {
    *violation += aj * solval;
    if aj > 0.0 && solval <= feastol {
        return;
    }
    if aj < 0.0 && solval >= upper - feastol {
        return;
    }
    *norm += aj * aj;
}

/// HiGHS: HighsCutGeneration.cpp::fast_floor (lines 502-502)
/// C++: `static double fast_floor(double x)`
pub fn fast_floor(x: f64) -> f64 {
    // C++ `(int64_t)x - (x < (int64_t)x)`: truncation toward zero, then minus one when the
    // truncation rounded up; the result is converted back to double.
    let t = x as i64;
    (t - ((x < t as f64) as i64)) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `5 x0 + 4 x1 + 3 x2 <= 8` over binaries with the cover `{x0, x1}` (weight 9, so
    /// `lambda = 1`), the state `determineCover` leaves. Hand-traced through the C++:
    /// abar = 4, S = [4, 8], x0 is in C+ and lifts to g(5) = 1, x1 gets 1, g(3) = 0, so the
    /// lifted cover is `x0 + x1 <= 1`.
    fn knapsack_cover_state() -> (HighsCutGeneration, Vec<f64>) {
        let cg = HighsCutGeneration {
            feastol: 1e-6,
            epsilon: 1e-9,
            rowlen: 3,
            rhs: HighsCDouble::from(8.0),
            lambda: HighsCDouble::from(1.0),
            cover: vec![1, 0],
            isintegral: vec![1, 1, 1],
            upper: vec![1.0, 1.0, 1.0],
            solval: vec![1.0, 1.0, 0.0],
            ..Default::default()
        };
        (cg, vec![5.0, 4.0, 3.0])
    }

    /// `determineCover(lpSol = true)` from scratch: x2 sits at 0 so it is skipped, x0 and x1
    /// sit at their upper bound and are partitioned in first; weight 9 > 8 gives lambda = 1.
    #[test]
    fn determine_cover_on_tiny_knapsack() {
        let (mut cg, mut vals) = knapsack_cover_state();
        cg.cover.clear();
        cg.lambda = HighsCDouble::from(0.0);
        let mut inds = vec![0, 1, 2];
        let ms = HighsMipSolver::default();
        assert!(cg.determine_cover(&ms, &mut vals, &mut inds, true));
        assert_eq!(cg.cover, vec![0, 1]);
        assert_eq!(f64::from(cg.lambda), 1.0);
        assert_eq!(f64::from(cg.coverweight), 9.0);
    }

    #[test]
    fn lifted_knapsack_cover_on_tiny_knapsack() {
        let (mut cg, mut vals) = knapsack_cover_state();
        cg.separate_lifted_knapsack_cover(&mut vals);
        assert_eq!(
            cg.cover,
            vec![0, 1],
            "cover sorted by decreasing coefficient"
        );
        assert_eq!(vals, vec![1.0, 1.0, 0.0]);
        assert_eq!(f64::from(cg.rhs), 1.0);
        assert!(cg.integral_support && cg.integral_coefficients);
    }

    /// Same knapsack plus a continuous `-y` (`y` complemented to its lower bound): the
    /// mixed-binary lifting keeps `-y`, caps the cover coefficients at lambda, and
    /// phi(3) = 0 for x2 (3 <= S[0] - lambda = 4): `x0 + x1 - y <= 1`.
    #[test]
    fn lifted_mixed_binary_cover_on_tiny_knapsack() {
        let (mut cg, mut vals) = knapsack_cover_state();
        cg.rowlen = 4;
        cg.isintegral.push(0);
        cg.upper.push(kHighsInf);
        cg.solval.push(0.0);
        vals.push(-1.0);
        assert!(cg.separate_lifted_mixed_binary_cover(&mut vals));
        assert_eq!(vals, vec![1.0, 1.0, 0.0, -1.0]);
        assert_eq!(f64::from(cg.rhs), 1.0);
        assert!(!cg.integral_support && !cg.integral_coefficients);
    }

    #[test]
    fn fast_floor_matches_cpp_truncation_formula() {
        for (x, want) in [
            (2.5, 2.0),
            (-2.5, -3.0),
            (-3.0, -3.0),
            (0.0, 0.0),
            (-0.5, -1.0),
        ] {
            assert_eq!(fast_floor(x), want, "fast_floor({x})");
        }
    }
}

#[cfg(test)]
mod h16_big_act_tests {
    use super::*;

    /// One integer column in [0, upper] with `vals = [val]`, feastol 1e-6, epsilon 1e-9.
    /// Returns (return value, has_unbounded_ints, has_general_ints).
    fn preprocess_int_case(on: bool, val: f64, upper: f64) -> (bool, bool, bool) {
        let mut ms = HighsMipSolver::default();
        let lp = &mut ms.orig_model_;
        lp.num_col_ = 1;
        lp.num_row_ = 0;
        lp.col_cost_ = vec![0.0];
        lp.col_lower_ = vec![0.0];
        lp.col_upper_ = vec![upper];
        lp.integrality_ = vec![HighsVarType::kInteger];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.h16_big_act = on;
        let mut lp_relaxation = HighsLpRelaxation::default();
        lp_relaxation.lpsolver.model_lp_.num_col_ = 1;
        let mut cg = HighsCutGeneration {
            feastol: 1e-6,
            epsilon: 1e-9,
            rowlen: 1,
            upper: vec![upper],
            solval: vec![0.0],
            ..Default::default()
        };
        let mut vals = vec![val];
        let mut inds = vec![0];
        let (mut has_unbounded_ints, mut has_general_ints, mut has_continuous) =
            (false, false, false);
        let ok = cg.preprocess_base_inequality(
            &ms,
            &lp_relaxation,
            &mut vals,
            &mut inds,
            &mut has_unbounded_ints,
            &mut has_general_ints,
            &mut has_continuous,
        );
        (ok, has_unbounded_ints, has_general_ints)
    }

    // HiGHS 1.16 #3267: HighsCutGeneration.cpp lines 963-965 (1.16.0): |vals * upper| > 2^52 -> reject.
    // vals = 0.5 is not rescaled by `scale`; 0.5 * 1e16 = 5e15 > 2^52.
    #[test]
    fn h16_big_act_on_preprocess_rejects_huge_product() {
        let (ok, unbounded, general) = preprocess_int_case(true, 0.5, 1e16);
        assert!(!ok);
        assert!(
            !unbounded && !general,
            "the rejection happens before the classification"
        );
    }

    #[test]
    fn h16_big_act_off_preprocess_keeps_huge_product() {
        // 1.15.1 has no product check: the integer is classified as a general integer.
        let (_, unbounded, general) = preprocess_int_case(false, 0.5, 1e16);
        assert!(!unbounded);
        assert!(general);
    }

    // HiGHS 1.16 #3267: HighsCutGeneration.cpp line 966 (1.16.0): `upper >= kHighsIInf` marks an
    // unbounded integer. upper = kHighsIInf is 2^31 - 1 here, so |vals * upper| < 2^52.
    #[test]
    fn h16_big_act_on_preprocess_unbounded_at_khighs_iinf() {
        let upper = kHighsIInf as f64;
        let (_, unbounded, general) = preprocess_int_case(true, 1e-3, upper);
        assert!(unbounded);
        assert!(general);
    }

    #[test]
    fn h16_big_act_off_preprocess_khighs_iinf_not_unbounded() {
        // 1.15.1 tests `upper == kHighsInf`, which a finite upper bound never satisfies.
        let upper = kHighsIInf as f64;
        let (_, unbounded, general) = preprocess_int_case(false, 1e-3, upper);
        assert!(!unbounded);
        assert!(general);
    }
}
