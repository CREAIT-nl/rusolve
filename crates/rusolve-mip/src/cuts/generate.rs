//! Single-row cut generation from a base inequality `Σ a_k y_k ≤ β` in the
//! bound-substituted space of [`super::transform`]: lifted knapsack cover,
//! lifted mixed-binary cover, lifted mixed-integer cover and complemented
//! mixed-integer rounding (c-MIR), the strongest-efficacy one wins.
//!
//! Sources:
//! * c-MIR: H. Marchand & L. A. Wolsey, "Aggregation and mixed integer
//!   rounding to solve MIPs", Operations Research 49(3):363–371 (2001) —
//!   delta trials, complementation flips, efficacy as the selection score.
//! * Lifted covers and superadditive lifting functions: Z. Gu, G. L.
//!   Nemhauser & M. W. P. Savelsbergh, "Lifted cover inequalities for 0-1
//!   integer programs: computation", INFORMS J. Computing 10(4):427–437
//!   (1998), and "Sequence independent lifting in mixed integer
//!   programming", J. Comb. Optim. 4:109–129 (2000); mixed-integer cover
//!   lifting per Letchford & Souli (2019) as implemented in HiGHS.
//! * The concrete algorithm, constants and safeguards follow HiGHS
//!   `HighsCutGeneration` (reference/HiGHS/highs/mip/HighsCutGeneration.cpp,
//!   MIT; read as a reference text and reimplemented): `determineCover`
//!   :25–129, `separateLiftedKnapsackCover` :131–214,
//!   `separateLiftedMixedBinaryCover` :216–278,
//!   `separateLiftedMixedIntegerCover` :280–500,
//!   `cmirCutGenerationHeuristic` :504–739 (f0 ∈ [0.005, 0.995], 1/(1−f0) ≤
//!   1e6, δ·2/4/8 trials), `scale` :741–758, `postprocessCut` :760–899,
//!   `preprocessBaseInequality` :901–1037 (row length cap 100 + 0.15·n),
//!   `generateCut` :1058–1197 (violation > 10·feastol),
//!   `tryGenerateCut` :1382–1482; `HighsDomain::tightenCoefficients`
//!   (HighsDomain.cpp:2668–2710); `HighsIntegers::integralScale/denominator`
//!   (util/HighsIntegers.h:87–203).
//!
//! Deviations: plain f64 instead of HiGHS' double-double accumulators; cover
//! ordering uses exact keys (Rust sorts need a total order) instead of
//! feastol-tolerant comparisons; the lifting loops carry iteration guards
//! that abandon the cut instead of looping on numerically odd input.

use super::relax::Relaxation;
use super::transform::TransformedLp;

/// HiGHS kHighsTiny.
const TINY: f64 = 1e-14;
const MAX_CMIR_SCALE: f64 = 1e6;
const F0_MIN: f64 = 0.005;
const F0_MAX: f64 = 0.995;
/// Bound on lifting-function inner iterations before the cut is dropped.
const LIFT_GUARD: u64 = 1_000_000;

/// Which inequality family produced a cut.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Method {
    KnapsackCover,
    MixedBinaryCover,
    MixedIntegerCover,
    Cmir,
}

pub(crate) struct GeneratedCut {
    pub inds: Vec<usize>,
    pub vals: Vec<f64>,
    pub rhs: f64,
    pub integral: bool,
    pub method: Method,
}

/// `x = m·2^e`, `m ∈ [0.5, 1)` — the exponent of C `frexp` (0 for 0/∞/NaN).
fn frexp_exp(x: f64) -> i32 {
    if x == 0.0 || !x.is_finite() {
        return 0;
    }
    let bits = x.abs().to_bits();
    let e = ((bits >> 52) & 0x7ff) as i32;
    if e == 0 {
        // subnormal
        return frexp_exp(x * exp2i(54)) - 54;
    }
    e - 1022
}

/// Exact 2^k for k in the normal range.
fn exp2i(k: i32) -> f64 {
    let k = k.clamp(-1022, 1023);
    f64::from_bits(((k + 1023) as u64) << 52)
}

fn fast_floor(x: f64) -> f64 {
    x.floor()
}

pub(crate) fn nearest_integer(x: f64) -> f64 {
    (x + 0.5f64.copysign(x)).trunc()
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let h = a % b;
        a = b;
        b = h;
    }
    a
}

/// Continued-fraction denominator of `x` up to `maxdenom`
/// (HighsIntegers::denominator).
fn denominator(x: f64, eps: f64, maxdenom: i64) -> i64 {
    let mut ai = x as i64;
    let mut m = [ai, 1i64, 1i64, 0i64];
    let mut xi = x;
    let mut fraction = xi - ai as f64;
    while fraction > eps {
        xi = 1.0 / fraction;
        if xi > (1u64 << 53) as f64 {
            break;
        }
        ai = xi as i64;
        let t = m[2] * ai + m[3];
        if t > maxdenom {
            break;
        }
        m[3] = m[2];
        m[2] = t;
        let t = m[0] * ai + m[1];
        m[1] = m[0];
        m[0] = t;
        fraction = xi - ai as f64;
    }
    ai = (maxdenom - m[3]) / m[2];
    m[1] += m[0] * ai;
    m[3] += m[2] * ai;
    let x0 = m[0] as f64 / m[2] as f64;
    let x1 = m[1] as f64 / m[3] as f64;
    let x = x.abs();
    if (x - x0).abs() < (x - x1).abs() {
        m[2]
    } else {
        m[3]
    }
}

/// Smallest scale making every value integral, or 0 (HighsIntegers::
/// integralScale).
pub(crate) fn integral_scale(vals: &[f64], deltadown: f64, deltaup: f64) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    let mut minval = vals[0];
    let mut maxval = vals[0];
    for &v in vals {
        if v.abs() < minval.abs() {
            minval = v;
        }
        if v.abs() > maxval.abs() {
            maxval = v;
        }
    }
    let mut expshift = 0;
    if minval < -deltadown || minval > deltaup {
        expshift = frexp_exp(minval);
    }
    expshift = (-expshift).max(0) + 3;
    let exp_max = frexp_exp(maxval).min(32);
    if exp_max + expshift > 32 {
        expshift = 32 - exp_max;
    }
    if !(0..=62).contains(&expshift) {
        return 0.0;
    }
    let mut denom: u64 = 75u64 << expshift;
    let mut startdenom: u64 = denom;
    let mut val = startdenom as f64 * vals[0];
    let mut downval = (val + deltaup).floor();
    let mut fraction = val - downval;
    if fraction > deltadown {
        denom = denom.saturating_mul(denominator(fraction, deltaup, 1000) as u64);
        val = denom as f64 * vals[0];
        downval = (val + deltaup).floor();
        fraction = val - downval;
        if fraction > deltadown {
            return 0.0;
        }
    }
    let mut currgcd = downval.abs() as u64;
    for &v in &vals[1..] {
        val = denom as f64 * v;
        downval = (val + deltaup).floor();
        fraction = val - downval;
        if fraction > deltadown {
            val = startdenom as f64 * v;
            fraction = val - val.floor();
            denom = denom.saturating_mul(denominator(fraction, deltaup, 1000) as u64);
            val = denom as f64 * v;
            downval = (val + deltaup).floor();
            fraction = val - downval;
            if fraction > deltadown {
                return 0.0;
            }
        }
        if currgcd != 1 {
            currgcd = gcd(currgcd, downval.abs() as u64);
            if denom > u32::MAX as u64 {
                if currgcd == 0 {
                    return 0.0;
                }
                denom /= currgcd;
                if startdenom != 1 {
                    startdenom /= gcd(currgcd, startdenom);
                }
                currgcd = 1;
            }
        }
    }
    if currgcd == 0 {
        return 0.0;
    }
    denom as f64 / currgcd as f64
}

fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Working state of one base inequality (HighsCutGeneration's members).
pub(crate) struct CutGen {
    feastol: f64,
    epsilon: f64,
    max_len: usize,
    rng: u64,
    inds: Vec<usize>,
    vals: Vec<f64>,
    upper: Vec<f64>,
    solval: Vec<f64>,
    isint: Vec<bool>,
    comp: Vec<bool>,
    rhs: f64,
    initial_scale: f64,
    cover: Vec<usize>,
    lambda: f64,
    integral_support: bool,
    integral_coefficients: bool,
    method: Method,
}

struct Flags {
    unbounded_ints: bool,
    general_ints: bool,
    continuous: bool,
}

impl CutGen {
    pub fn new(num_col: usize, feastol: f64, epsilon: f64, seed: u64) -> Self {
        Self {
            feastol,
            epsilon,
            max_len: 100 + (0.15 * num_col as f64) as usize,
            rng: seed,
            inds: Vec::new(),
            vals: Vec::new(),
            upper: Vec::new(),
            solval: Vec::new(),
            isint: Vec::new(),
            comp: Vec::new(),
            rhs: 0.0,
            initial_scale: 1.0,
            cover: Vec::new(),
            lambda: 0.0,
            integral_support: false,
            integral_coefficients: false,
            method: Method::Cmir,
        }
    }

    /// Generates a cut from `Σ vals·x[inds] ≤ rhs` (inds over columns and
    /// slacks) violated by more than 10·feastol at `relax.x`, in original
    /// column space. `only_initial_scale` restricts c-MIR to the row's own
    /// scale (no δ candidates from coefficients, no ×2/4/8 trials), as for
    /// mod-k aggregations whose scale is already the right one
    /// (HighsCutGeneration.cpp:536, :563, :614).
    pub fn generate(
        &mut self,
        relax: &Relaxation,
        trans: &mut TransformedLp,
        inds: Vec<usize>,
        vals: Vec<f64>,
        rhs: f64,
        only_initial_scale: bool,
    ) -> Option<GeneratedCut> {
        let t = trans.transform(relax, inds, vals, rhs, true, false)?;
        let n = t.inds.len();
        self.inds = t.inds;
        self.vals = t.vals;
        self.upper = t.upper;
        self.solval = t.solval;
        self.rhs = t.rhs;
        self.comp = vec![false; n];
        let flags = self.preprocess(relax)?;
        if !flags.unbounded_ints && !t.integers_positive {
            for i in 0..self.inds.len() {
                if self.vals[i] > 0.0 || !self.isint[i] {
                    continue;
                }
                self.flip(i);
            }
        }
        if !self.try_generate(&flags, 10.0 * self.feastol, only_initial_scale) {
            return None;
        }
        for i in 0..self.inds.len() {
            if self.comp[i] {
                self.flip(i);
            }
        }
        let (mut ci, mut cv) = (Vec::new(), Vec::new());
        for (&k, &v) in self.inds.iter().zip(&self.vals) {
            if v != 0.0 {
                ci.push(k);
                cv.push(v);
            }
        }
        let (inds, vals, rhs) = trans.untransform(relax, &ci, &cv, self.rhs, false)?;
        self.inds = inds;
        self.vals = vals;
        self.rhs = rhs;
        if !self.postprocess(relax) {
            return None;
        }
        let viol: f64 = self
            .inds
            .iter()
            .zip(&self.vals)
            .map(|(&j, &v)| v * relax.x[j])
            .sum::<f64>()
            - self.rhs;
        if viol <= 10.0 * self.feastol {
            return None;
        }
        self.tighten_coefficients(relax);
        Some(GeneratedCut {
            inds: std::mem::take(&mut self.inds),
            vals: std::mem::take(&mut self.vals),
            rhs: self.rhs,
            integral: self.integral_support && self.integral_coefficients,
            method: self.method,
        })
    }

    fn flip(&mut self, i: usize) {
        debug_assert!(self.upper[i].is_finite());
        self.comp[i] = !self.comp[i];
        self.solval[i] = self.upper[i] - self.solval[i];
        self.rhs -= self.upper[i] * self.vals[i];
        self.vals[i] = -self.vals[i];
    }

    fn update_viol_norm(&self, i: usize, aj: f64, viol: &mut f64, norm: &mut f64) {
        *viol += aj * self.solval[i];
        if aj > 0.0 && self.solval[i] <= self.feastol {
            return;
        }
        if aj < 0.0 && self.solval[i] >= self.upper[i] - self.feastol {
            return;
        }
        *norm += aj * aj;
    }

    /// Scales by the power of two that brings `val` to about 1 (at most
    /// ×1024); returns the factor.
    fn scale(&mut self, val: f64) -> f64 {
        let expshift = (-frexp_exp(val)).min(10);
        let f = exp2i(expshift);
        self.rhs *= f;
        for v in &mut self.vals {
            *v *= f;
        }
        f
    }

    fn retain_nonzero(&mut self) {
        let mut k = 0;
        for i in 0..self.vals.len() {
            if self.vals[i] != 0.0 {
                self.inds.swap(k, i);
                self.vals.swap(k, i);
                self.upper.swap(k, i);
                self.solval.swap(k, i);
                self.isint.swap(k, i);
                self.comp.swap(k, i);
                k += 1;
            }
        }
        self.inds.truncate(k);
        self.vals.truncate(k);
        self.upper.truncate(k);
        self.solval.truncate(k);
        self.isint.truncate(k);
        self.comp.truncate(k);
    }

    fn preprocess(&mut self, relax: &Relaxation) -> Option<Flags> {
        let feastol = self.feastol;
        let mut flags = Flags {
            unbounded_ints: false,
            general_ints: false,
            continuous: false,
        };
        let mut maxact = -feastol;
        let max_abs = self.vals.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        self.initial_scale = self.scale(max_abs);
        let n = self.inds.len();
        self.isint = vec![false; n];
        let mut nonzeros = n;
        for i in 0..n {
            self.isint[i] = relax.is_integral(self.inds[i]) && self.vals[i].abs() > 10.0 * feastol;
            if !self.isint[i] {
                if self.upper[i] < 2.0 * self.solval[i] {
                    self.flip(i);
                }
                if self.vals[i] > 0.0 || self.vals[i].abs() * self.upper[i] <= 10.0 * feastol {
                    if self.vals[i] < 0.0 {
                        if self.upper[i] == f64::INFINITY {
                            return None;
                        }
                        self.rhs -= self.vals[i] * self.upper[i];
                    }
                    if self.vals[i] != 0.0 {
                        nonzeros -= 1;
                    }
                    self.vals[i] = 0.0;
                    continue;
                }
                flags.continuous = true;
            } else {
                if self.upper[i] == f64::INFINITY {
                    flags.unbounded_ints = true;
                    flags.general_ints = true;
                } else if self.upper[i] != 1.0 {
                    flags.general_ints = true;
                }
                if self.vals[i] > 0.0 {
                    maxact += self.vals[i] * self.upper[i];
                }
            }
        }
        if nonzeros > self.max_len {
            let num_cancel = nonzeros - self.max_len;
            let mut cancel: Vec<usize> = (0..n)
                .filter(|&i| {
                    self.vals[i] != 0.0 && {
                        let slack = if self.vals[i] > 0.0 {
                            self.solval[i]
                        } else {
                            self.upper[i] - self.solval[i]
                        };
                        slack <= feastol
                    }
                })
                .collect();
            if cancel.len() < num_cancel {
                return None;
            }
            cancel.sort_unstable_by(|&a, &b| {
                self.vals[a]
                    .abs()
                    .total_cmp(&self.vals[b].abs())
                    .then(a.cmp(&b))
            });
            for &j in &cancel[..num_cancel] {
                if self.vals[j] < 0.0 {
                    self.rhs -= self.vals[j] * self.upper[j];
                } else {
                    maxact -= self.vals[j] * self.upper[j];
                }
                self.vals[j] = 0.0;
            }
        }
        self.retain_nonzero();
        if maxact > self.rhs { Some(flags) } else { None }
    }

    fn determine_cover(&mut self) -> bool {
        let feastol = self.feastol;
        if self.rhs <= 10.0 * feastol {
            return false;
        }
        let n = self.inds.len();
        let (mut at_upper, mut rest): (Vec<usize>, Vec<usize>) = (0..n)
            .filter(|&j| self.isint[j] && self.solval[j] > feastol)
            .partition(|&j| self.solval[j] >= self.upper[j] - feastol);
        let mut coverweight: f64 = at_upper.iter().map(|&j| self.vals[j] * self.upper[j]).sum();
        self.rng = splitmix(self.rng);
        let r = self.rng;
        rest.sort_unstable_by(|&i, &j| {
            let bi = self.upper[i] > 1.5;
            let bj = self.upper[j] > 1.5;
            bi.cmp(&bj)
                .then((self.solval[j] * self.vals[j]).total_cmp(&(self.solval[i] * self.vals[i])))
                .then(self.vals[j].total_cmp(&self.vals[i]))
                .then(splitmix(self.inds[j] as u64 ^ r).cmp(&splitmix(self.inds[i] as u64 ^ r)))
        });
        let minlambda = (10.0 * feastol).max(feastol * self.rhs.abs());
        let mut k = 0;
        while k < rest.len() {
            if coverweight - self.rhs > minlambda {
                break;
            }
            let j = rest[k];
            coverweight += self.vals[j] * self.upper[j];
            k += 1;
        }
        at_upper.extend_from_slice(&rest[..k]);
        if at_upper.is_empty() {
            return false;
        }
        self.lambda = coverweight - self.rhs;
        if self.lambda <= minlambda {
            return false;
        }
        self.cover = at_upper;
        true
    }

    fn sort_cover_desc(&mut self) {
        let vals = &self.vals;
        self.cover
            .sort_unstable_by(|&a, &b| vals[b].total_cmp(&vals[a]).then(a.cmp(&b)));
    }

    // Parallel index arrays (cover / S / coverflag) read clearer indexed.
    #[allow(clippy::needless_range_loop)]
    fn lifted_knapsack_cover(&mut self) {
        let feastol = self.feastol;
        let coversize = self.cover.len();
        let n = self.inds.len();
        self.sort_cover_desc();
        let cover = &self.cover;
        let vals = &mut self.vals;
        let mut abar = vals[cover[0]];
        let mut sigma = self.lambda;
        for i in 1..coversize {
            let delta = abar - vals[cover[i]];
            let kdelta = i as f64 * delta;
            if kdelta < sigma {
                abar = vals[cover[i]];
                sigma -= kdelta;
            } else {
                abar -= sigma * (1.0 / i as f64);
                sigma = 0.0;
                break;
            }
        }
        if sigma > 0.0 {
            abar = self.rhs / coversize as f64;
        }
        let mut s = vec![0.0; coversize];
        let mut coverflag = vec![0i8; n];
        let mut sum = 0.0;
        let mut cplussize = 0usize;
        for i in 0..coversize {
            sum += abar.min(vals[cover[i]]);
            s[i] = sum;
            if vals[cover[i]] > abar + feastol {
                cplussize += 1;
                coverflag[cover[i]] = 1;
            } else {
                coverflag[cover[i]] = -1;
            }
        }
        let mut halfintegral = false;
        let eps = self.epsilon;
        let mut g = |z: f64| -> f64 {
            let hfrac = z / abar;
            let mut coef = 0.0;
            let mut h = (hfrac + 0.5).floor() as i64;
            if h != 0 && (hfrac - h as f64).abs() * abar.max(1.0) <= eps && h < cplussize as i64 {
                halfintegral = true;
                coef = 0.5;
            }
            h = (h - 1).max(0);
            while h < coversize as i64 && z > s[h as usize] + feastol {
                h += 1;
            }
            coef + h as f64
        };
        self.rhs = coversize as f64 - 1.0;
        for i in 0..n {
            if vals[i] == 0.0 {
                continue;
            }
            vals[i] = if coverflag[i] == -1 { 1.0 } else { g(vals[i]) };
        }
        if halfintegral {
            self.rhs *= 2.0;
            for v in vals.iter_mut() {
                *v *= 2.0;
            }
        }
        self.integral_support = true;
        self.integral_coefficients = true;
        self.method = Method::KnapsackCover;
    }

    #[allow(clippy::needless_range_loop)]
    fn lifted_mixed_binary_cover(&mut self) -> bool {
        self.integral_support = false;
        self.integral_coefficients = false;
        let coversize = self.cover.len();
        if coversize == 0 {
            return false;
        }
        let n = self.inds.len();
        let mut coverflag = vec![false; n];
        for &c in &self.cover {
            coverflag[c] = true;
        }
        self.sort_cover_desc();
        let lambda = self.lambda;
        let mut s = vec![0.0; coversize];
        let mut sum = 0.0;
        let mut p = coversize;
        for i in 0..coversize {
            let v = self.vals[self.cover[i]];
            if v - lambda <= self.epsilon {
                p = i;
                break;
            }
            sum += v;
            s[i] = sum;
        }
        if p == 0 {
            return false;
        }
        let phi = |a: f64| -> f64 {
            for i in 0..p {
                if a <= s[i] - lambda {
                    return i as f64 * lambda;
                }
                if a <= s[i] {
                    return (i + 1) as f64 * lambda + (a - s[i]);
                }
            }
            p as f64 * lambda + (a - s[p - 1])
        };
        self.rhs = -lambda;
        self.integral_support = true;
        for i in 0..n {
            if !self.isint[i] {
                if self.vals[i] < 0.0 {
                    self.integral_support = false;
                } else {
                    self.vals[i] = 0.0;
                }
                continue;
            }
            if coverflag[i] {
                self.vals[i] = self.vals[i].min(lambda);
                self.rhs += self.vals[i];
            } else {
                self.vals[i] = phi(self.vals[i]);
            }
        }
        self.method = Method::MixedBinaryCover;
        true
    }

    fn lifted_mixed_integer_cover(&mut self) -> bool {
        let feastol = self.feastol;
        self.integral_support = false;
        self.integral_coefficients = false;
        let n = self.inds.len();
        let mut coverflag = vec![false; n];
        for &c in &self.cover {
            coverflag[c] = true;
        }
        self.sort_cover_desc();
        let mut cover = self.cover.clone();
        let coversize = cover.len();
        let vals = &self.vals;
        let upper = &self.upper;
        let mut a: Vec<f64> = Vec::with_capacity(coversize);
        let mut u: Vec<f64> = Vec::with_capacity(coversize + 1);
        let mut m: Vec<f64> = Vec::with_capacity(coversize + 1);
        let (mut usum, mut msum) = (0.0, 0.0);
        for &i in &cover {
            u.push(usum);
            m.push(msum);
            a.push(vals[i]);
            usum += upper[i];
            msum += upper[i] * vals[i];
        }
        u.push(usum);
        m.push(msum);

        let lambda = self.lambda;
        let mut lpos: Option<usize> = None;
        let mut best_cplusend = 0usize;
        let mut best_val = 0.0;
        let mut best_at_upper = true;
        for i in 0..coversize {
            let j = cover[i];
            let ub = upper[j];
            let at_upper = self.solval[j] >= ub - feastol;
            if at_upper && !best_at_upper {
                continue;
            }
            let mju = ub * vals[j];
            let mu = mju - lambda;
            if mu <= 10.0 * feastol || vals[j].abs() < 1000.0 * feastol {
                continue;
            }
            let mudival = mu / vals[j];
            let fr = (mudival - mudival.trunc()).abs();
            if fr.min(1.0 - fr) <= feastol {
                continue;
            }
            let eta = mudival.ceil();
            let cplusthreshold = (ub - eta + 1.0) * vals[j];
            let cplusend = cover.partition_point(|&c| vals[c] >= cplusthreshold);
            let mut mcplus = m[cplusend];
            if i < cplusend {
                mcplus -= mju;
            }
            let jl_val = mcplus + eta * vals[j];
            if jl_val > best_val || (!at_upper && best_at_upper) {
                lpos = Some(i);
                best_cplusend = cplusend;
                best_val = jl_val;
                best_at_upper = at_upper;
            }
        }
        let Some(lpos) = lpos else {
            return false;
        };
        let l = cover[lpos];
        let al = vals[l];
        let upperl = upper[l];
        let mlu = upperl * al;
        let mu = mlu - lambda;
        a.truncate(best_cplusend);
        cover.truncate(best_cplusend);
        u.truncate(best_cplusend + 1);
        m.truncate(best_cplusend + 1);
        if lpos < best_cplusend {
            a.remove(lpos);
            cover.remove(lpos);
            u.remove(lpos + 1);
            m.remove(lpos + 1);
            for i in (lpos + 1)..best_cplusend {
                u[i] -= upperl;
                m[i] -= mlu;
            }
        }
        let cplussize = a.len();
        let mudival = mu / al;
        let eta = mudival.ceil();
        let r = (mu - mudival.floor() * al).max(0.0);
        let ulminusetaplusone = upperl - eta + 1.0;
        let cplusthreshold = ulminusetaplusone * al;
        let kmin = (eta - upperl - 0.5).floor() as i64;
        let kmax = (upperl - eta + 0.5).floor() as i64;
        let mut guard = 0u64;

        let phi_l = |x: f64, guard: &mut u64| -> Option<f64> {
            let mut k = ((x / al) as i64).min(-1);
            while k >= kmin {
                *guard += 1;
                if *guard > LIFT_GUARD {
                    return None;
                }
                let kf = k as f64;
                if x >= kf * al + r {
                    return Some(x - (kf + 1.0) * r);
                }
                if x >= kf * al {
                    return Some(kf * (al - r));
                }
                k -= 1;
            }
            Some(kmin as f64 * (al - r))
        };
        let gamma_l = |z: f64, guard: &mut u64| -> Option<f64> {
            for i in 0..cplussize {
                let upperi = upper[cover[i]] as i64;
                for h in 0..=upperi {
                    let hf = h as f64;
                    let mih = m[i] + hf * a[i];
                    let uih = u[i] + hf;
                    let mihplusdeltai = mih + a[i] - cplusthreshold;
                    *guard += 1;
                    if *guard > LIFT_GUARD {
                        return None;
                    }
                    if z <= mihplusdeltai {
                        return Some(uih * ulminusetaplusone * (al - r));
                    }
                    let mut k = ((z - mihplusdeltai) / al) as i64 - 1;
                    while k <= kmax {
                        *guard += 1;
                        if *guard > LIFT_GUARD {
                            return None;
                        }
                        let kf = k as f64;
                        if z <= mihplusdeltai + kf * al + r {
                            return Some((uih * ulminusetaplusone + kf) * (al - r));
                        }
                        if z <= mihplusdeltai + (kf + 1.0) * al {
                            return Some(
                                uih * ulminusetaplusone * (al - r) + z - mih - a[i]
                                    + cplusthreshold
                                    - (kf + 1.0) * r,
                            );
                        }
                        k += 1;
                    }
                }
            }
            let mut p = ((z - m[cplussize]) / al) as i64 - 1;
            loop {
                *guard += 1;
                if *guard > LIFT_GUARD {
                    return None;
                }
                let pf = p as f64;
                if z <= m[cplussize] + pf * al + r {
                    return Some((u[cplussize] * ulminusetaplusone + pf) * (al - r));
                }
                if z <= m[cplussize] + (pf + 1.0) * al {
                    return Some(
                        u[cplussize] * ulminusetaplusone * (al - r) + z
                            - m[cplussize]
                            - (pf + 1.0) * r,
                    );
                }
                p += 1;
            }
        };

        let mut new_vals = self.vals.clone();
        let mut rhs = (upperl - eta) * r - lambda;
        let mut support = true;
        for i in 0..n {
            if new_vals[i] == 0.0 {
                continue;
            }
            if !self.isint[i] {
                if new_vals[i] < 0.0 {
                    support = false;
                } else {
                    new_vals[i] = 0.0;
                }
                continue;
            }
            if coverflag[i] {
                let Some(v) = phi_l(-new_vals[i], &mut guard) else {
                    return false;
                };
                new_vals[i] = -v;
                rhs += new_vals[i] * upper[i];
            } else {
                let Some(v) = gamma_l(new_vals[i], &mut guard) else {
                    return false;
                };
                new_vals[i] = v;
            }
        }
        self.vals = new_vals;
        self.rhs = rhs;
        self.integral_support = support;
        self.method = Method::MixedIntegerCover;
        true
    }

    /// Efficacy of the c-MIR cut for step `delta` (or `None` if f0 is out of
    /// range).
    fn cmir_efficacy(
        &self,
        delta: f64,
        integerinds: &[usize],
        cont_contrib: f64,
        cont_sqnorm: f64,
    ) -> Option<f64> {
        let scale = 1.0 / delta;
        let scalrhs = self.rhs * scale;
        let downrhs = fast_floor(scalrhs);
        let f0 = scalrhs - downrhs;
        if !(F0_MIN..=F0_MAX).contains(&f0) {
            return None;
        }
        let k = 1.0 / (1.0 - f0);
        if k > MAX_CMIR_SCALE {
            return None;
        }
        let contscale = scale * k;
        let mut sqrnorm = contscale * contscale * cont_sqnorm;
        let mut viol = contscale * cont_contrib - downrhs;
        for &j in integerinds {
            let scalaj = self.vals[j] * scale;
            let downaj = fast_floor(scalaj + TINY);
            let fj = scalaj - downaj;
            let aj = downaj + ((fj - f0) * k).max(0.0);
            self.update_viol_norm(j, aj, &mut viol, &mut sqrnorm);
        }
        Some(viol / sqrnorm.sqrt())
    }

    fn cmir(&mut self, min_efficacy: f64, only_initial_scale: bool) -> bool {
        let feastol = self.feastol;
        self.integral_support = false;
        self.integral_coefficients = false;
        let n = self.inds.len();
        let (mut cont_contrib, mut cont_sqnorm) = (0.0, 0.0);
        let mut deltas: Vec<f64> = Vec::with_capacity(n + 3);
        let mut integerinds: Vec<usize> = Vec::with_capacity(n);
        let mut maxabsdelta = 0.0f64;
        for i in 0..n {
            if self.isint[i] {
                integerinds.push(i);
                if self.upper[i] < 2.0 * self.solval[i] {
                    self.flip(i);
                }
                if only_initial_scale {
                    continue;
                }
                if self.solval[i] > feastol {
                    let delta = self.vals[i].abs();
                    if delta <= 1e-4 || delta == maxabsdelta {
                        continue;
                    }
                    maxabsdelta = maxabsdelta.max(delta);
                    deltas.push(delta);
                }
            } else {
                let v = self.vals[i];
                self.update_viol_norm(i, v, &mut cont_contrib, &mut cont_sqnorm);
            }
        }
        if cont_sqnorm == 0.0 && deltas.len() > 1 {
            let int_scale = integral_scale(&deltas, feastol, TINY);
            if int_scale != 0.0 && int_scale <= 1e4 {
                let scalrhs = self.rhs * int_scale;
                let f0 = scalrhs - fast_floor(scalrhs);
                if (F0_MIN..=F0_MAX).contains(&f0) {
                    deltas.push(1.0 / int_scale);
                }
            }
        }
        deltas.push(self.initial_scale.min(1.0));
        if !only_initial_scale {
            deltas.push(maxabsdelta + self.initial_scale.min(1.0));
        }
        deltas.sort_unstable_by(f64::total_cmp);
        let mut cur = deltas[0];
        for d in deltas.iter_mut().skip(1) {
            if *d - cur <= 10.0 * feastol {
                *d = 0.0;
            } else {
                cur = *d;
            }
        }
        deltas.retain(|&d| d != 0.0);
        let mut bestdelta = -1.0;
        let mut besteff = min_efficacy;
        for &delta in &deltas {
            if let Some(eff) = self.cmir_efficacy(delta, &integerinds, cont_contrib, cont_sqnorm)
                && eff > besteff
            {
                bestdelta = delta;
                besteff = eff;
            }
        }
        if bestdelta == -1.0 {
            return false;
        }
        for k in (1..=3).filter(|_| !only_initial_scale) {
            let delta = bestdelta * f64::from(1u32 << k);
            if let Some(eff) = self.cmir_efficacy(delta, &integerinds, cont_contrib, cont_sqnorm)
                && eff > besteff
            {
                bestdelta = delta;
                besteff = eff;
            }
        }
        for &k in &integerinds {
            if self.upper[k] == f64::INFINITY || self.solval[k] <= feastol {
                continue;
            }
            self.flip(k);
            match self.cmir_efficacy(bestdelta, &integerinds, cont_contrib, cont_sqnorm) {
                Some(eff) if eff > besteff => besteff = eff,
                _ => self.flip(k),
            }
        }
        let scale = 1.0 / bestdelta;
        let scalrhs = self.rhs * scale;
        let downrhs = scalrhs.floor();
        let f0 = scalrhs - downrhs;
        let k = 1.0 / (1.0 - f0);
        self.rhs = downrhs * bestdelta;
        self.integral_support = true;
        for j in 0..n {
            if self.vals[j] == 0.0 {
                continue;
            }
            if !self.isint[j] {
                if self.vals[j] > 0.0 {
                    self.vals[j] = 0.0;
                } else {
                    self.vals[j] *= k;
                    self.integral_support = false;
                }
            } else {
                let scalaj = scale * self.vals[j];
                let downaj = (scalaj + TINY).floor();
                let fj = scalaj - downaj;
                let mut aj = downaj;
                if fj > f0 {
                    aj += (fj - f0) * k;
                }
                self.vals[j] = aj * bestdelta;
            }
        }
        self.method = Method::Cmir;
        true
    }

    fn try_generate(&mut self, flags: &Flags, min_efficacy: f64, only_initial_scale: bool) -> bool {
        if flags.unbounded_ints {
            return self.cmir(min_efficacy, only_initial_scale);
        }
        let tmp_vals = self.vals.clone();
        let tmp_comp = self.comp.clone();
        let tmp_solval = self.solval.clone();
        let mut tmp_rhs = self.rhs;
        let mut success = false;
        if self.determine_cover() {
            if !flags.continuous && !flags.general_ints {
                self.lifted_knapsack_cover();
                success = true;
            } else if flags.general_ints {
                success = self.lifted_mixed_integer_cover();
            } else {
                success = self.lifted_mixed_binary_cover();
            }
        }
        let mut min_mir = min_efficacy;
        let (mut save_support, mut save_coefs, mut save_method) = (false, false, Method::Cmir);
        if success {
            save_support = self.integral_support;
            save_coefs = self.integral_coefficients;
            save_method = self.method;
            let (mut viol, mut sqrnorm) = (-self.rhs, 0.0);
            for i in 0..self.inds.len() {
                let v = self.vals[i];
                self.update_viol_norm(i, v, &mut viol, &mut sqrnorm);
            }
            let eff = viol / sqrnorm.sqrt();
            if eff <= min_efficacy || eff.is_nan() {
                success = false;
                self.rhs = tmp_rhs;
            } else {
                min_mir += eff;
                std::mem::swap(&mut tmp_rhs, &mut self.rhs);
            }
        } else {
            self.rhs = tmp_rhs;
        }
        let lifted_vals = std::mem::replace(&mut self.vals, tmp_vals);
        if self.cmir(min_mir, only_initial_scale) {
            return true;
        }
        if success {
            self.rhs = tmp_rhs;
            self.comp = tmp_comp;
            self.solval = tmp_solval;
            self.vals = lifted_vals;
            self.integral_support = save_support;
            self.integral_coefficients = save_coefs;
            self.method = save_method;
            return true;
        }
        false
    }

    fn postprocess(&mut self, relax: &Relaxation) -> bool {
        let feastol = self.feastol;
        if self.rhs < 0.0 && self.rhs > -self.epsilon {
            self.rhs = 0.0;
        }
        if self.integral_support && self.integral_coefficients {
            self.drop_zero_entries();
            return true;
        }
        let mut max_abs = self.vals.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let min_coef = 100.0 * feastol * max_abs.max(1e-3);
        self.integral_support = true;
        for i in 0..self.inds.len() {
            let v = self.vals[i];
            if v == 0.0 {
                continue;
            }
            let j = self.inds[i];
            if v.abs() <= min_coef {
                if v < 0.0 {
                    if relax.upper[j] == f64::INFINITY {
                        return false;
                    }
                    self.rhs -= relax.upper[j] * v;
                } else {
                    if relax.lower[j] == f64::NEG_INFINITY {
                        return false;
                    }
                    self.rhs -= relax.lower[j] * v;
                }
                self.vals[i] = 0.0;
                continue;
            }
            if !relax.integer[j] {
                self.integral_support = false;
            }
        }
        self.drop_zero_entries();
        if self.inds.is_empty() {
            return false;
        }
        if self.integral_support {
            let intscale = integral_scale(&self.vals, feastol, self.epsilon);
            let mut scale_smallest_to_one = true;
            if intscale != 0.0 && intscale * max_abs.max(1.0) <= (1u64 << 52) as f64 {
                self.rhs *= intscale;
                max_abs = nearest_integer(max_abs * intscale);
                for i in 0..self.inds.len() {
                    let scaleval = intscale * self.vals[i];
                    let intval = nearest_integer(scaleval);
                    let delta = scaleval - intval;
                    self.vals[i] = intval;
                    let j = self.inds[i];
                    if delta < 0.0 {
                        if relax.upper[j] == f64::INFINITY {
                            return false;
                        }
                        self.rhs -= delta * relax.upper[j];
                    } else {
                        if relax.lower[j] == f64::NEG_INFINITY {
                            return false;
                        }
                        self.rhs -= delta * relax.lower[j];
                    }
                }
                self.rhs = (self.rhs + feastol).floor();
                if intscale * max_abs * feastol < 0.5 {
                    scale_smallest_to_one = false;
                    self.integral_coefficients = true;
                }
            }
            // Integer rounding can zero a coefficient.
            self.drop_zero_entries();
            if self.inds.is_empty() {
                return false;
            }
            if scale_smallest_to_one {
                let min_abs = self.vals.iter().fold(f64::INFINITY, |m, v| m.min(v.abs()));
                self.scale(min_abs - self.epsilon);
            }
        } else {
            self.scale(max_abs - self.epsilon);
        }
        true
    }

    fn drop_zero_entries(&mut self) {
        let mut k = 0;
        for i in 0..self.vals.len() {
            if self.vals[i] != 0.0 {
                self.inds.swap(k, i);
                self.vals.swap(k, i);
                k += 1;
            }
        }
        self.inds.truncate(k);
        self.vals.truncate(k);
    }

    fn tighten_coefficients(&mut self, relax: &Relaxation) {
        let mut maxact = 0.0;
        for (&j, &v) in self.inds.iter().zip(&self.vals) {
            if v > 0.0 {
                if relax.upper[j] == f64::INFINITY {
                    return;
                }
                maxact += relax.upper[j] * v;
            } else {
                if relax.lower[j] == f64::NEG_INFINITY {
                    return;
                }
                maxact += relax.lower[j] * v;
            }
        }
        let maxabscoef = maxact - self.rhs;
        if maxabscoef > self.feastol {
            let mut upper = self.rhs;
            let mut tightened = false;
            for i in 0..self.inds.len() {
                let j = self.inds[i];
                if !relax.integer[j] {
                    continue;
                }
                if self.vals[i] > maxabscoef {
                    upper -= (self.vals[i] - maxabscoef) * relax.upper[j];
                    self.vals[i] = maxabscoef;
                    tightened = true;
                } else if self.vals[i] < -maxabscoef {
                    upper += (-self.vals[i] - maxabscoef) * relax.lower[j];
                    self.vals[i] = -maxabscoef;
                    tightened = true;
                }
            }
            if tightened {
                self.rhs = upper;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frexp_and_exp2_match_c() {
        assert_eq!(frexp_exp(1.0), 1);
        assert_eq!(frexp_exp(0.5), 0);
        assert_eq!(frexp_exp(0.75), 0);
        assert_eq!(frexp_exp(3.0), 2);
        assert_eq!(frexp_exp(-3.0), 2);
        assert_eq!(frexp_exp(1e-310), -1029);
        assert_eq!(exp2i(10), 1024.0);
        assert_eq!(exp2i(-3), 0.125);
    }

    #[test]
    fn integral_scale_finds_common_denominator() {
        let s = integral_scale(&[0.5, 1.0 / 3.0, 2.0], 1e-6, 1e-14);
        assert!(s > 0.0);
        for v in [0.5, 1.0 / 3.0, 2.0] {
            let x = v * s;
            assert!((x - x.round()).abs() < 1e-6, "{v} * {s}");
        }
        assert!((s - 6.0).abs() < 1e-9, "{s}");
        assert_eq!(
            integral_scale(&[std::f64::consts::PI, 1.0], 1e-9, 1e-14),
            0.0
        );
    }

    #[test]
    fn denominator_of_simple_fractions() {
        assert_eq!(denominator(1.0 / 3.0, 1e-9, 1000), 3);
        assert_eq!(denominator(0.75, 1e-9, 1000), 4);
    }
}
