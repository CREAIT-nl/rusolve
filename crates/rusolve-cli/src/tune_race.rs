//! Statistics of `rusolve tune`: the racing test, a bootstrap interval and a
//! small deterministic RNG. Pure functions, no I/O.
//!
//! ## Racing: what is followed, and where it was read
//!
//! * F-Race: Birattari, Stuetzle, Paquete, Varrentrapp, "A Racing Algorithm for Configuring
//!   Metaheuristics", GECCO 2002. Candidates are evaluated block by block (one block is one
//!   (instance, seed) pair, run by every live candidate); after each block a Friedman test on
//!   the per-block ranks decides whether any candidate is worse than the best, and the
//!   candidates that are statistically worse are dropped for good.
//! * irace: Lopez-Ibanez, Dubois-Lacoste, Perez Caceres, Birattari, Stuetzle, "The irace
//!   package: Iterated racing for automatic algorithm configuration", Operations Research
//!   Perspectives 3 (2016) 43-58. Defaults used here are the ones the irace reference page
//!   (mlopez-ibanez.github.io/irace/reference/irace.html, read for this change) prints for a
//!   scenario: `testType = "friedman"`, `confidence = 0.95`, `firstTest = 5` (no test before
//!   five blocks), `eachTest = 1` (test after every block), `elitist = TRUE` (the best
//!   candidates found so far are kept; here: the baseline is never dropped during the race, so
//!   the final comparison against it always has data). Only the racing test and these defaults
//!   live in this file; the iterated sampling of new candidates from the elite of a round (the
//!   `nbIterations` / `elitist` half of irace) is in `tune.rs::tune`, and it is our own
//!   simplification of irace's model- and test-instances machinery, not a port of it.
//! * The test and the post-hoc comparison are Conover's: W. J. Conover, Practical
//!   Nonparametric Statistics, 3rd ed., Wiley 1999, section 5.8 (Friedman test, statistic T2
//!   compared with the F distribution, and the pairwise rank-sum comparison with a Student t
//!   quantile). The formulas are written out at [`friedman`]. The page above does not print
//!   the post-hoc formula, so these come from the textbook, not from irace's source.
//! * The F and t tail probabilities use the regularised incomplete beta function by the
//!   continued fraction of Numerical Recipes (Press et al., 3rd ed., section 6.4, `betacf`,
//!   modified Lentz) and the Lanczos log-gamma (NR section 6.1). Their values are tested
//!   against closed forms and tabulated quantiles below.
//! * Bootstrap: Efron, "Bootstrap methods: another look at the jackknife", Annals of
//!   Statistics 7 (1979), percentile interval.

/// Log-gamma, Lanczos approximation (g = 7, n = 9; Numerical Recipes section 6.1).
pub fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = G[0];
    let t = x + 7.5;
    for (i, g) in G.iter().enumerate().skip(1) {
        a += g / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..300 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-14 {
            break;
        }
    }
    h
}

/// Regularised incomplete beta function I_x(a, b).
pub fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

/// P(F > f) for F with (d1, d2) degrees of freedom.
pub fn f_sf(f: f64, d1: f64, d2: f64) -> f64 {
    if f <= 0.0 {
        return 1.0;
    }
    betai(d2 / 2.0, d1 / 2.0, d2 / (d2 + d1 * f))
}

/// Two-sided tail P(|T| > t) for Student t with `df` degrees of freedom.
pub fn t_sf2(t: f64, df: f64) -> f64 {
    betai(df / 2.0, 0.5, df / (df + t * t))
}

/// The t quantile with two-sided tail probability `alpha` (found by bisection).
pub fn t_crit(df: f64, alpha: f64) -> f64 {
    let (mut lo, mut hi) = (0.0, 1e4);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if t_sf2(mid, df) > alpha {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Ranks of one block, 1 = lowest cost; ties get the average rank.
pub fn ranks(row: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..row.len()).collect();
    order.sort_by(|&a, &b| row[a].total_cmp(&row[b]));
    let mut r = vec![0.0; row.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j + 1 < order.len() && row[order[j + 1]] == row[order[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0 + 1.0;
        for &k in &order[i..=j] {
            r[k] = avg;
        }
        i = j + 1;
    }
    r
}

/// Mean rank of every column of a blocks x candidates table (1 = cheapest in a block, ties
/// averaged). The elite of an iterated-racing round is ranked by this, which is what irace's
/// reports rank its elite configurations by as well. Empty table, empty result.
pub fn mean_ranks(table: &[Vec<f64>]) -> Vec<f64> {
    let Some(first) = table.first() else {
        return vec![];
    };
    let mut sums = vec![0.0f64; first.len()];
    for row in table {
        for (j, r) in ranks(row).into_iter().enumerate() {
            sums[j] += r;
        }
    }
    sums.iter().map(|s| s / table.len() as f64).collect()
}

/// Outcome of a Friedman test on a blocks x candidates table.
#[derive(Debug, Clone)]
pub struct Friedman {
    pub p_value: f64,
    pub rank_sums: Vec<f64>,
    /// Rank-sum difference that the post-hoc comparison calls significant.
    pub critical_difference: f64,
}

/// Conover's Friedman test (Practical Nonparametric Statistics, 3rd ed., 5.8).
///
/// With `b` blocks, `k` candidates, rank sums `R_j`, `A1 = sum of squared ranks` and
/// `C1 = b k (k+1)^2 / 4`:
///
/// * `T1 = (k-1) (sum_j R_j^2 - b C1) / (A1 - C1)`, and `T2 = (b-1) T1 / (b (k-1) - T1)`,
///   which follows F(k-1, (b-1)(k-1)) under the hypothesis that all candidates are alike;
/// * candidates i and j differ at level alpha when
///   `|R_i - R_j| > t(1 - alpha/2, (b-1)(k-1)) * sqrt(2 (b A1 - sum_j R_j^2) / ((b-1)(k-1)))`.
///
/// `None` when the test is undefined: fewer than two blocks or candidates, or every block a
/// complete tie.
pub fn friedman(table: &[Vec<f64>], alpha: f64) -> Option<Friedman> {
    let b = table.len();
    let k = table.first()?.len();
    if b < 2 || k < 2 {
        return None;
    }
    let (bf, kf) = (b as f64, k as f64);
    let mut rank_sums = vec![0.0; k];
    let mut a1 = 0.0;
    for row in table {
        for (j, r) in ranks(row).into_iter().enumerate() {
            rank_sums[j] += r;
            a1 += r * r;
        }
    }
    let c1 = bf * kf * (kf + 1.0).powi(2) / 4.0;
    let sum_r2: f64 = rank_sums.iter().map(|r| r * r).sum();
    if (a1 - c1).abs() < 1e-12 {
        return None;
    }
    let t1 = (kf - 1.0) * (sum_r2 - bf * c1) / (a1 - c1);
    let denom = bf * (kf - 1.0) - t1;
    let p_value = if denom <= 1e-12 {
        0.0
    } else {
        f_sf((bf - 1.0) * t1 / denom, kf - 1.0, (bf - 1.0) * (kf - 1.0))
    };
    let df = (bf - 1.0) * (kf - 1.0);
    let spread = (2.0 * (bf * a1 - sum_r2) / df).max(0.0).sqrt();
    Some(Friedman {
        p_value,
        rank_sums,
        critical_difference: t_crit(df, alpha) * spread,
    })
}

/// Columns of `table` that are significantly worse than the best column after a significant
/// Friedman test at level `alpha` (F-Race elimination rule). Empty when the test is not
/// significant or undefined.
pub fn race_drops(table: &[Vec<f64>], alpha: f64) -> Vec<usize> {
    let Some(f) = friedman(table, alpha) else {
        return vec![];
    };
    if f.p_value >= alpha {
        return vec![];
    }
    let best = f.rank_sums.iter().cloned().fold(f64::INFINITY, f64::min);
    (0..f.rank_sums.len())
        .filter(|&j| f.rank_sums[j] - best > f.critical_difference)
        .collect()
}

/// splitmix64: small, seedable and reproducible across platforms.
pub struct Rng(pub u64);

impl Rng {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    /// Fisher-Yates.
    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            v.swap(i, self.below(i + 1));
        }
    }
}

/// Percentile bootstrap of the mean of `xs` (log-ratios): (mean, lower, upper) at 95%.
pub fn bootstrap_mean_ci(xs: &[f64], resamples: usize, seed: u64) -> (f64, f64, f64) {
    let n = xs.len();
    let mean = xs.iter().sum::<f64>() / n.max(1) as f64;
    if n < 2 {
        return (mean, mean, mean);
    }
    let mut rng = Rng(seed);
    let mut means: Vec<f64> = (0..resamples)
        .map(|_| (0..n).map(|_| xs[rng.below(n)]).sum::<f64>() / n as f64)
        .collect();
    means.sort_by(f64::total_cmp);
    let at = |q: f64| means[((q * (resamples - 1) as f64).round() as usize).min(resamples - 1)];
    (mean, at(0.025), at(0.975))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn special_functions_match_known_values() {
        // Gamma(5) = 24, Gamma(1/2) = sqrt(pi)
        assert!((ln_gamma(5.0) - 24f64.ln()).abs() < 1e-12);
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-12);
        // I_x(1, 1) = x; I_x(2, 1) = x^2
        assert!((betai(1.0, 1.0, 0.3) - 0.3).abs() < 1e-12);
        assert!((betai(2.0, 1.0, 0.3) - 0.09).abs() < 1e-12);
        // two-sided 5% t quantiles from the standard table
        for (df, q) in [(1.0, 12.706), (5.0, 2.571), (10.0, 2.228), (30.0, 2.042)] {
            assert!((t_crit(df, 0.05) - q).abs() < 2e-3, "df {df}");
        }
        // 5% critical value of F(3, 20) is 3.098
        assert!((f_sf(3.098, 3.0, 20.0) - 0.05).abs() < 1e-3);
    }

    #[test]
    fn ranks_average_ties() {
        assert_eq!(ranks(&[3.0, 1.0, 1.0, 2.0]), vec![4.0, 1.5, 1.5, 3.0]);
    }

    #[test]
    fn mean_ranks_average_the_per_block_ranks() {
        // column 1 is cheapest in both blocks, column 0 most expensive in both
        let table = vec![vec![3.0, 1.0, 2.0], vec![3.0, 1.0, 2.0]];
        assert_eq!(mean_ranks(&table), vec![3.0, 1.0, 2.0]);
        // a tie in one block is averaged: block 1 ranks (1.5,1.5), block 2 ranks (2,1)
        let ties = vec![vec![1.0, 1.0], vec![2.0, 1.0]];
        assert_eq!(mean_ranks(&ties), vec![1.75, 1.25]);
        assert!(mean_ranks(&[]).is_empty());
        assert!(mean_ranks(&[]).is_empty());
    }

    #[test]
    fn dominated_candidate_is_dropped_and_equals_are_kept() {
        // column 2 is always worst by a wide margin; columns 0 and 1 trade places
        let table: Vec<Vec<f64>> = (0..8)
            .map(|i| {
                let e = if i % 2 == 0 { 0.1 } else { -0.1 };
                vec![10.0 + e, 10.0 - e, 50.0 + i as f64]
            })
            .collect();
        assert_eq!(race_drops(&table, 0.05), vec![2]);
        // pure noise: nothing dropped
        let noise: Vec<Vec<f64>> = (0..8)
            .map(|i| {
                let mut r = Rng(i as u64 + 1);
                (0..4).map(|_| r.below(1000) as f64).collect()
            })
            .collect();
        assert!(race_drops(&noise, 0.05).is_empty());
        // all ties: the test is undefined, nothing dropped
        assert!(race_drops(&vec![vec![1.0, 1.0, 1.0]; 6], 0.05).is_empty());
    }

    #[test]
    fn bootstrap_interval_brackets_the_mean() {
        let xs: Vec<f64> = (0..20).map(|i| -0.3 + 0.01 * (i % 5) as f64).collect();
        let (m, lo, hi) = bootstrap_mean_ci(&xs, 500, 7);
        assert!(lo <= m && m <= hi && hi < 0.0);
    }
}
