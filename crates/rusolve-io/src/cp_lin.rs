//! MIP linearization of the [`CpModel`] objects: the engine-free fallback of the CP layer.
//!
//! `linearize` derives binaries and rows for a solve from the user's objects; the objects stay in
//! the model and the checker evaluates them, never the derived rows. New columns are appended
//! after the model's `n` columns and new rows after its rows, so cutting a solution back to the
//! user's `n` columns is exact (the same contract as [`crate::gc_lin`], whose [`Extension`] it
//! returns).
//!
//! Where each formulation comes from:
//! - Interval `start + size = end`, optional presence `p`: the equality as two big-M rows
//!   relaxed by `p`, with `M` read off the column bounds (Williams, "Model Building in
//!   Mathematical Programming", 5th ed., ch. 9; the indicator-to-big-M step is the one of Bonami,
//!   Lodi, Tramontani, Wiese, Math. Prog. 151, 2015).
//! - No-overlap: one disjunctive binary per pair, `end_i <= start_j` or `end_j <= start_i`, big-M
//!   from the bounds, relaxed by the presence columns. Manne, "On the job-shop scheduling
//!   problem", Oper. Res. 8 (1960); Applegate and Cook, "A computational study of the job-shop
//!   scheduling problem", ORSA J. Comput. 3 (1991). A pair that the bounds already order is
//!   dropped.
//! - Cumulative: the time-indexed formulation, one binary per (interval, start time) and one
//!   capacity row per time point. Pritsker, Watters and Wolfe, "Multiproject scheduling with
//!   limited resources: a zero-one programming approach", Management Science 16 (1969). It needs
//!   constant sizes and demands and a horizon of at most [`CpLinOptions::horizon_max`] time
//!   points; otherwise a typed refusal says the CP-SAT engine handles it.
//! - All-different: the direct (value-binary) encoding, one binary per (variable, value), each
//!   variable takes exactly one, each value at most once. Williams, ch. 9 (assignment
//!   structure); Walsh, "SAT v CSP", CP 2000 (the direct encoding).
//! - Element `target = values[index]`: one selection binary per admissible index, `index` and
//!   `target` as the selection-weighted sums (big-M rows when a value is a column). Williams,
//!   ch. 9; Hooker, "Integrated Methods for Optimization", 2nd ed. (2012), element constraint.
//! - Table: one binary per allowed tuple, each column the tuple-weighted sum (the same
//!   references). A negated table forbids each tuple with a disjunction over the columns
//!   (`x_j <= v - 1` or `x_j >= v + 1`, one binary per side, big-M from the bounds).
//!
//! A valid model the formulation cannot take is a typed refusal ([`CpErrorKind::Unsupported`]):
//! a big-M above `bigm_max`, a cumulative horizon above the limit, a variable size or demand in a
//! cumulative, a no-overlap size column that may be 0, a domain too large to enumerate.

use crate::cp::{CpError, CpModel, CpTerm};
use crate::gc_lin::{Extension, NewCol, NewRow, Route};
use crate::model::VariableBound;

/// Options of the linearization.
#[derive(Debug, Clone, Copy)]
pub struct CpLinOptions {
    /// Largest big-M a derived row may carry.
    pub bigm_max: f64,
    /// Largest horizon (time points) of one cumulative in the time-indexed formulation.
    pub horizon_max: usize,
    /// Largest number of enumerated binaries per object (value, tuple or start-time binaries).
    pub binaries_max: usize,
}

impl Default for CpLinOptions {
    fn default() -> Self {
        Self {
            bigm_max: 1e6,
            horizon_max: 1000,
            binaries_max: 100_000,
        }
    }
}

/// What the linearization reads of the user's model.
pub struct CpLinInput<'a> {
    pub col_lower: &'a [f64],
    pub col_upper: &'a [f64],
    /// Whether each column is integer.
    pub col_integer: &'a [bool],
    pub col_names: &'a [String],
}

struct Lin<'a> {
    n: usize,
    lo: &'a [f64],
    up: &'a [f64],
    opt: CpLinOptions,
    ext: Extension,
}

fn ceil_i(v: f64) -> i64 {
    (v - 1e-9).ceil() as i64
}

fn floor_i(v: f64) -> i64 {
    (v + 1e-9).floor() as i64
}

impl Lin<'_> {
    fn lb(&self, j: usize) -> i64 {
        ceil_i(self.lo[j])
    }

    fn ub(&self, j: usize) -> i64 {
        floor_i(self.up[j])
    }

    fn tb(&self, t: &CpTerm) -> (i64, i64) {
        match *t {
            CpTerm::Const(c) => (c, c),
            CpTerm::Col(j) => (self.lb(j), self.ub(j)),
        }
    }

    fn bin(&mut self, name: String) -> usize {
        self.ext.cols.push(NewCol {
            name,
            lower: 0.0,
            upper: 1.0,
            integer: true,
        });
        self.n + self.ext.cols.len() - 1
    }

    /// `lower <= sum terms <= upper`; repeated columns are summed, zeros dropped.
    fn row(&mut self, name: String, terms: Vec<(usize, f64)>, lower: f64, upper: f64) {
        let mut merged: Vec<(usize, f64)> = Vec::with_capacity(terms.len());
        for (j, c) in terms {
            match merged.iter_mut().find(|(k, _)| *k == j) {
                Some(e) => e.1 += c,
                None => merged.push((j, c)),
            }
        }
        merged.retain(|&(_, c)| c != 0.0);
        self.ext.rows.push(NewRow {
            name,
            terms: merged,
            lower,
            upper,
        });
    }

    /// A row no point satisfies (an object whose bounds already contradict it).
    fn infeasible(&mut self, name: String) {
        self.row(name, Vec::new(), 1.0, 1.0);
    }

    fn route(&mut self, object: &str, route: &'static str, max_big_m: f64) {
        self.ext.routes.push(Route {
            object: object.to_string(),
            route,
            max_big_m,
        });
    }

    /// A big-M above the limit is a typed refusal, never a silent huge coefficient.
    fn big(&self, obj: &str, m: i64) -> Result<f64, CpError> {
        let m = m as f64;
        if m > self.opt.bigm_max {
            return Err(CpError::unsupported(
                obj,
                format!(
                    "the linearization needs a big-M of {m} (from the variable bounds), above the limit {}; tighten the bounds, or use the CP-SAT engine",
                    self.opt.bigm_max
                ),
            ));
        }
        Ok(m)
    }

    fn check_count(&self, obj: &str, what: &str, count: usize) -> Result<(), CpError> {
        if count > self.opt.binaries_max {
            return Err(CpError::unsupported(
                obj,
                format!(
                    "the linearization needs {count} {what} binaries, above the limit {}; use the CP-SAT engine",
                    self.opt.binaries_max
                ),
            ));
        }
        Ok(())
    }
}

/// Derives the MIP rows and binaries of `cp`. `Err` with [`CpErrorKind::Invalid`] for a model
/// that fails [`CpModel::validate`], with [`CpErrorKind::Unsupported`] for one the fallback cannot
/// take.
pub fn linearize(
    cp: &CpModel,
    inp: &CpLinInput<'_>,
    opt: CpLinOptions,
) -> Result<Extension, CpError> {
    let bounds: Vec<VariableBound> = (0..inp.col_lower.len())
        .map(|j| VariableBound {
            lower: inp.col_lower[j],
            upper: inp.col_upper[j],
            integer: inp.col_integer.get(j).copied().unwrap_or(false),
        })
        .collect();
    cp.validate(&bounds, inp.col_names)?;
    let mut l = Lin {
        n: inp.col_lower.len(),
        lo: inp.col_lower,
        up: inp.col_upper,
        opt,
        ext: Extension::default(),
    };
    for it in &cp.intervals {
        l.interval(it)?;
    }
    for o in &cp.no_overlap {
        l.no_overlap(cp, &o.name, &o.intervals)?;
    }
    for o in &cp.cumulative {
        l.cumulative(cp, o)?;
    }
    for o in &cp.all_different {
        l.all_different(&o.name, &o.cols)?;
    }
    for o in &cp.element {
        l.element(o)?;
    }
    for o in &cp.table {
        l.table(o)?;
    }
    Ok(l.ext)
}

impl Lin<'_> {
    fn interval(&mut self, it: &crate::cp::Interval) -> Result<(), CpError> {
        let obj = format!("interval '{}'", it.name);
        let (szlo, szhi) = self.tb(&it.size);
        let mut vars = vec![(it.start, 1.0), (it.end, -1.0)];
        let mut c = 0i64;
        match it.size {
            CpTerm::Const(v) => c = v,
            CpTerm::Col(j) => vars.push((j, 1.0)),
        }
        let name = format!("iv_{}", it.name);
        let Some(p) = it.presence else {
            self.row(name, vars, -(c as f64), -(c as f64));
            self.route(&obj, "interval_equality", 0.0);
            return Ok(());
        };
        // E = start + size - end ranges over [lo_e, up_e]; E = 0 is enforced when p = 1.
        let lo_e = self.lb(it.start) + szlo - self.ub(it.end);
        let up_e = self.ub(it.start) + szhi - self.lb(it.end);
        let mut m = 0;
        if up_e > 0 {
            // E <= up_e (1 - p)
            let big = self.big(&obj, up_e)?;
            let mut t = vars.clone();
            t.push((p, big));
            self.row(format!("{name}_le"), t, f64::NEG_INFINITY, big - c as f64);
            m = m.max(up_e);
        }
        if lo_e < 0 {
            // E >= lo_e (1 - p)
            let big = self.big(&obj, -lo_e)?;
            let mut t = vars;
            t.push((p, -big));
            self.row(format!("{name}_ge"), t, -big - c as f64, f64::INFINITY);
            m = m.max(-lo_e);
        }
        self.route(&obj, "interval_big_m", m as f64);
        Ok(())
    }

    fn no_overlap(&mut self, cp: &CpModel, name: &str, list: &[usize]) -> Result<(), CpError> {
        let obj = format!("no_overlap '{name}'");
        // a constant zero size overlaps nothing; the checker skips it as well
        let live: Vec<usize> = list
            .iter()
            .copied()
            .filter(|&i| cp.intervals[i].size != CpTerm::Const(0))
            .collect();
        for &i in &live {
            if let CpTerm::Col(j) = cp.intervals[i].size
                && self.lb(j) < 1
            {
                return Err(CpError::unsupported(
                    &obj,
                    format!(
                        "interval '{}' has a size variable that can be 0; a zero-size interval overlaps nothing, but the disjunctive rows would still order it. Give the size a lower bound of at least 1, or use the CP-SAT engine",
                        cp.intervals[i].name
                    ),
                ));
            }
        }
        let mut max_m = 0;
        let mut pairs = 0usize;
        for (a, &i) in live.iter().enumerate() {
            for &j in &live[a + 1..] {
                let (p, q) = (&cp.intervals[i], &cp.intervals[j]);
                // i before j needs end_i - start_j <= m1; j before i needs end_j - start_i <= m2
                let m1 = self.ub(p.end) - self.lb(q.start);
                let m2 = self.ub(q.end) - self.lb(p.start);
                if m1 <= 0 || m2 <= 0 {
                    continue; // the bounds order the pair
                }
                let (b1, b2) = (self.big(&obj, m1)?, self.big(&obj, m2)?);
                max_m = max_m.max(m1).max(m2);
                let y = self.bin(format!("no_{name}_{}_{}", p.name, q.name));
                pairs += 1;
                let pres: Vec<usize> = [p.presence, q.presence].into_iter().flatten().collect();
                let np = pres.len() as f64;
                // y = 1: i before j.  end_i - start_j + m1 y + m1 sum(p) <= m1 (1 + np)
                let mut t = vec![(p.end, 1.0), (q.start, -1.0), (y, b1)];
                t.extend(pres.iter().map(|&k| (k, b1)));
                self.row(
                    format!("no_{name}_{}_{}_a", p.name, q.name),
                    t,
                    f64::NEG_INFINITY,
                    b1 * (1.0 + np),
                );
                // y = 0: j before i.  end_j - start_i - m2 y + m2 sum(p) <= m2 np
                let mut t = vec![(q.end, 1.0), (p.start, -1.0), (y, -b2)];
                t.extend(pres.iter().map(|&k| (k, b2)));
                self.row(
                    format!("no_{name}_{}_{}_b", p.name, q.name),
                    t,
                    f64::NEG_INFINITY,
                    b2 * np,
                );
            }
        }
        self.route(
            &obj,
            if pairs == 0 {
                "no_overlap_redundant"
            } else {
                "no_overlap_big_m"
            },
            max_m as f64,
        );
        Ok(())
    }

    fn cumulative(&mut self, cp: &CpModel, o: &crate::cp::Cumulative) -> Result<(), CpError> {
        let obj = format!("cumulative '{}'", o.name);
        struct Act<'x> {
            name: &'x str,
            start: usize,
            pres: Option<usize>,
            size: i64,
            demand: i64,
            a: i64,
            b: i64,
        }
        let fixed = |me: &Self, t: &CpTerm, what: &str, who: &str| -> Result<i64, CpError> {
            let (lo, hi) = me.tb(t);
            if lo == hi {
                Ok(lo)
            } else {
                Err(CpError::unsupported(
                    &obj,
                    format!(
                        "{who} has a variable {what} [{lo}, {hi}]; the time-indexed formulation needs a constant {what}. The CP-SAT engine takes variable sizes and demands"
                    ),
                ))
            }
        };
        let mut acts: Vec<Act> = Vec::new();
        for (k, &i) in o.intervals.iter().enumerate() {
            let it = &cp.intervals[i];
            let who = format!("interval '{}'", it.name);
            let size = fixed(self, &it.size, "size", &who)?;
            let demand = fixed(self, &o.demands[k], "demand", &who)?;
            if size == 0 || demand == 0 {
                continue;
            }
            let a = self.lb(it.start).max(self.lb(it.end) - size);
            let b = self.ub(it.start).min(self.ub(it.end) - size);
            if a > b {
                // the interval cannot be present within its bounds
                match it.presence {
                    Some(p) => self.row(
                        format!("cum_{}_{}_absent", o.name, it.name),
                        vec![(p, 1.0)],
                        0.0,
                        0.0,
                    ),
                    None => self.infeasible(format!("cum_{}_{}_window", o.name, it.name)),
                }
                continue;
            }
            acts.push(Act {
                name: &it.name,
                start: it.start,
                pres: it.presence,
                size,
                demand,
                a,
                b,
            });
        }
        if acts.is_empty() {
            self.route(&obj, "cumulative_redundant", 0.0);
            return Ok(());
        }
        let tmin = acts.iter().map(|x| x.a).min().unwrap();
        let tmax = acts.iter().map(|x| x.b + x.size).max().unwrap();
        let horizon = (tmax - tmin) as usize;
        if horizon > self.opt.horizon_max {
            return Err(CpError::unsupported(
                &obj,
                format!(
                    "the time-indexed formulation needs a horizon of {horizon} time points (from the start and end bounds), above the limit {}; tighten the bounds, or use the CP-SAT engine",
                    self.opt.horizon_max
                ),
            ));
        }
        let total: usize = acts.iter().map(|x| (x.b - x.a + 1) as usize).sum();
        self.check_count(&obj, "start-time", total)?;
        // x[k][tau - a]: interval k starts at tau
        let mut x: Vec<Vec<usize>> = Vec::new();
        for ac in &acts {
            let mut v = Vec::new();
            for tau in ac.a..=ac.b {
                v.push(self.bin(format!("cum_{}_{}_{tau}", o.name, ac.name)));
            }
            // sum x = presence (1 when always present)
            let mut t: Vec<(usize, f64)> = v.iter().map(|&c| (c, 1.0)).collect();
            match ac.pres {
                Some(p) => {
                    t.push((p, -1.0));
                    self.row(format!("cum_{}_{}_one", o.name, ac.name), t, 0.0, 0.0);
                }
                None => self.row(format!("cum_{}_{}_one", o.name, ac.name), t, 1.0, 1.0),
            }
            // start = sum tau x when present, free within its bounds otherwise
            let mut t: Vec<(usize, f64)> = vec![(ac.start, 1.0)];
            t.extend(v.iter().zip(ac.a..).map(|(&c, tau)| (c, -(tau as f64))));
            match ac.pres {
                None => self.row(format!("cum_{}_{}_start", o.name, ac.name), t, 0.0, 0.0),
                Some(p) => {
                    let (sl, su) = (self.lb(ac.start), self.ub(ac.start));
                    let bu = self.big(&obj, su.abs())?;
                    let bl = self.big(&obj, sl.abs())?;
                    let mut t1 = t.clone();
                    t1.push((p, bu.copysign(su as f64)));
                    // start - sum + su p <= su
                    self.row(
                        format!("cum_{}_{}_start_le", o.name, ac.name),
                        t1,
                        f64::NEG_INFINITY,
                        su as f64,
                    );
                    let mut t2 = t;
                    t2.push((p, bl.copysign(sl as f64)));
                    // start - sum + sl p >= sl
                    self.row(
                        format!("cum_{}_{}_start_ge", o.name, ac.name),
                        t2,
                        sl as f64,
                        f64::INFINITY,
                    );
                }
            }
            x.push(v);
        }
        let (cap_lo, cap_term) = match o.capacity {
            CpTerm::Const(c) => (c, None),
            CpTerm::Col(j) => (self.lb(j), Some(j)),
        };
        for t in tmin..tmax {
            let mut terms: Vec<(usize, f64)> = Vec::new();
            let mut worst = 0i64;
            for (k, ac) in acts.iter().enumerate() {
                let lo = ac.a.max(t - ac.size + 1);
                let hi = ac.b.min(t);
                if lo > hi {
                    continue;
                }
                worst += ac.demand;
                for tau in lo..=hi {
                    terms.push((x[k][(tau - ac.a) as usize], ac.demand as f64));
                }
            }
            if terms.is_empty() || (cap_term.is_none() && worst <= cap_lo) {
                continue;
            }
            let upper = match cap_term {
                Some(j) => {
                    terms.push((j, -1.0));
                    0.0
                }
                None => cap_lo as f64,
            };
            self.row(
                format!("cum_{}_t{t}", o.name),
                terms,
                f64::NEG_INFINITY,
                upper,
            );
        }
        self.route(&obj, "cumulative_time_indexed", 0.0);
        Ok(())
    }

    fn all_different(&mut self, name: &str, cols: &[usize]) -> Result<(), CpError> {
        let obj = format!("all_different '{name}'");
        let mut distinct: Vec<usize> = cols.to_vec();
        distinct.sort_unstable();
        distinct.dedup();
        let total: i64 = distinct
            .iter()
            .map(|&j| (self.ub(j) - self.lb(j) + 1).max(0))
            .sum();
        self.check_count(&obj, "value", total as usize)?;
        // value -> [(position in `cols`, binary)]
        let mut by_value: std::collections::BTreeMap<i64, Vec<usize>> = Default::default();
        for &j in &distinct {
            let (lo, hi) = (self.lb(j), self.ub(j));
            if lo > hi {
                self.infeasible(format!("ad_{name}_{j}_empty"));
                continue;
            }
            let mut sel = vec![(j, 1.0)];
            let mut one = Vec::new();
            for v in lo..=hi {
                let z = self.bin(format!("ad_{name}_{j}_{v}"));
                sel.push((z, -(v as f64)));
                one.push((z, 1.0));
                // a column listed twice takes the binary twice: 2 z <= 1 forbids it
                for _ in 0..cols.iter().filter(|&&c| c == j).count() {
                    by_value.entry(v).or_default().push(z);
                }
            }
            self.row(format!("ad_{name}_{j}_one"), one, 1.0, 1.0);
            self.row(format!("ad_{name}_{j}_val"), sel, 0.0, 0.0);
        }
        for (v, zs) in by_value {
            if zs.len() > 1 {
                let t = zs.into_iter().map(|z| (z, 1.0)).collect();
                self.row(format!("ad_{name}_v{v}"), t, f64::NEG_INFINITY, 1.0);
            }
        }
        self.route(&obj, "all_different_value_binaries", 0.0);
        Ok(())
    }

    fn element(&mut self, o: &crate::cp::Element) -> Result<(), CpError> {
        let obj = format!("element '{}'", o.name);
        let n = o.values.len() as i64;
        let (lo, hi) = (self.lb(o.index).max(0), self.ub(o.index).min(n - 1));
        if lo > hi {
            self.infeasible(format!("el_{}_range", o.name));
            self.route(&obj, "element_selection", 0.0);
            return Ok(());
        }
        self.check_count(&obj, "selection", (hi - lo + 1) as usize)?;
        let w: Vec<(i64, usize)> = (lo..=hi)
            .map(|k| (k, self.bin(format!("el_{}_{k}", o.name))))
            .collect();
        self.row(
            format!("el_{}_one", o.name),
            w.iter().map(|&(_, c)| (c, 1.0)).collect(),
            1.0,
            1.0,
        );
        let mut t = vec![(o.index, 1.0)];
        t.extend(w.iter().map(|&(k, c)| (c, -(k as f64))));
        self.row(format!("el_{}_index", o.name), t, 0.0, 0.0);
        let mut max_m = 0;
        if o.values.iter().all(|v| matches!(v, CpTerm::Const(_))) {
            // target = sum values[k] w_k
            let mut t = vec![(o.target, 1.0)];
            for &(k, c) in &w {
                if let CpTerm::Const(v) = o.values[k as usize] {
                    t.push((c, -(v as f64)));
                }
            }
            self.row(format!("el_{}_target", o.name), t, 0.0, 0.0);
        } else {
            let (tlo, thi) = (self.lb(o.target), self.ub(o.target));
            for &(k, c) in &w {
                let v = o.values[k as usize];
                let (vlo, vhi) = self.tb(&v);
                // D = target - values[k] must be 0 when w_k = 1
                let (up_d, lo_d) = (thi - vlo, tlo - vhi);
                let (mut vars, off) = (vec![(o.target, 1.0)], 0i64);
                let off = match v {
                    CpTerm::Const(cv) => off - cv,
                    CpTerm::Col(j) => {
                        vars.push((j, -1.0));
                        off
                    }
                };
                if up_d > 0 {
                    let big = self.big(&obj, up_d)?;
                    max_m = max_m.max(up_d);
                    let mut t = vars.clone();
                    t.push((c, big));
                    self.row(
                        format!("el_{}_{k}_le", o.name),
                        t,
                        f64::NEG_INFINITY,
                        big - off as f64,
                    );
                }
                if lo_d < 0 {
                    let big = self.big(&obj, -lo_d)?;
                    max_m = max_m.max(-lo_d);
                    let mut t = vars.clone();
                    t.push((c, -big));
                    self.row(
                        format!("el_{}_{k}_ge", o.name),
                        t,
                        -big - off as f64,
                        f64::INFINITY,
                    );
                }
            }
        }
        self.route(&obj, "element_selection", max_m as f64);
        Ok(())
    }

    fn table(&mut self, o: &crate::cp::Table) -> Result<(), CpError> {
        let obj = format!("table '{}'", o.name);
        let doms: Vec<(i64, i64)> = o.cols.iter().map(|&j| (self.lb(j), self.ub(j))).collect();
        // a tuple outside the column box can never be taken
        let inside = |t: &Vec<i64>| {
            t.iter()
                .zip(&doms)
                .all(|(&v, &(lo, hi))| lo <= v && v <= hi)
        };
        let tuples: Vec<&Vec<i64>> = o.tuples.iter().filter(|t| inside(t)).collect();
        self.check_count(
            &obj,
            "tuple",
            tuples.len() * if o.negated { o.cols.len() * 2 } else { 1 },
        )?;
        if !o.negated {
            if tuples.is_empty() {
                self.infeasible(format!("tb_{}_empty", o.name));
                self.route(&obj, "table_tuple_binaries", 0.0);
                return Ok(());
            }
            let t: Vec<usize> = (0..tuples.len())
                .map(|k| self.bin(format!("tb_{}_{k}", o.name)))
                .collect();
            self.row(
                format!("tb_{}_one", o.name),
                t.iter().map(|&c| (c, 1.0)).collect(),
                1.0,
                1.0,
            );
            for (pos, &j) in o.cols.iter().enumerate() {
                let mut r = vec![(j, 1.0)];
                for (k, tup) in tuples.iter().enumerate() {
                    r.push((t[k], -(tup[pos] as f64)));
                }
                self.row(format!("tb_{}_c{pos}", o.name), r, 0.0, 0.0);
            }
            self.route(&obj, "table_tuple_binaries", 0.0);
            return Ok(());
        }
        // negated: every forbidden tuple differs in at least one column
        let mut max_m = 0;
        for (k, tup) in tuples.iter().enumerate() {
            let mut any: Vec<(usize, f64)> = Vec::new();
            for (pos, &j) in o.cols.iter().enumerate() {
                let (lo, hi) = doms[pos];
                let v = tup[pos];
                if lo < v {
                    // a = 1 -> x_j <= v - 1:  x_j + (hi - v + 1) a <= hi
                    let m = hi - v + 1;
                    let big = self.big(&obj, m)?;
                    max_m = max_m.max(m);
                    let a = self.bin(format!("tn_{}_{k}_{pos}_lt", o.name));
                    self.row(
                        format!("tn_{}_{k}_{pos}_lt", o.name),
                        vec![(j, 1.0), (a, big)],
                        f64::NEG_INFINITY,
                        hi as f64,
                    );
                    any.push((a, 1.0));
                }
                if v < hi {
                    // b = 1 -> x_j >= v + 1:  x_j - (v + 1 - lo) b >= lo
                    let m = v + 1 - lo;
                    let big = self.big(&obj, m)?;
                    max_m = max_m.max(m);
                    let b = self.bin(format!("tn_{}_{k}_{pos}_gt", o.name));
                    self.row(
                        format!("tn_{}_{k}_{pos}_gt", o.name),
                        vec![(j, 1.0), (b, -big)],
                        lo as f64,
                        f64::INFINITY,
                    );
                    any.push((b, 1.0));
                }
            }
            // no side available: the box is the forbidden tuple itself, `any` empty >= 1 fails
            self.row(format!("tn_{}_{k}_any", o.name), any, 1.0, f64::INFINITY);
        }
        self.route(&obj, "table_negated_disjunctive", max_m as f64);
        Ok(())
    }
}
