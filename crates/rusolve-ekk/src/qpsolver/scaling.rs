//! HiGHS source: qpsolver/scaling.hpp, qpsolver/scaling.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! NOTE: `largestpoweroftwo` returns `powf(1.0, l)` in HiGHS
//! (always 1.0); keep `1.0f32.powf(...)` semantics exactly as written in C++
//! (`powf` is the float overload: compute in `f32`, widen to `f64`).

use crate::prelude::*;
use crate::qpsolver::*;

/// HiGHS: scaling.cpp::largestpoweroftwo (lines 13-17)
fn largestpoweroftwo(value: f64) -> f64 {
    let l: f64 = value.log2();
    let il: HighsInt = l as HighsInt;
    // PORT NOTE: C++ powf(1.0, il) is always 1.0 (float overload, 1.0^anything = 1.0)
    (1.0f32).powf(il as f32) as f64
}

/// HiGHS: scaling.cpp::scale_rows (lines 19-58)
fn scale_rows(rt: &mut Runtime) {
    if !rt.settings.rowscaling {
        return;
    }

    // std::map operator[] -> `*map.entry(key).or_insert(0.0)` (inserts a
    // value-initialised 0.0 for a missing key, as C++ does).
    let mut maxabscoefperrow: BTreeMap<HighsInt, f64> = BTreeMap::new();
    let mut row: HighsInt = 0;
    while row < rt.scaled.num_con {
        *maxabscoefperrow.entry(row).or_insert(0.0) = 0.0;
        row += 1;
    }

    let mut var: HighsInt = 0;
    while var < rt.scaled.num_var {
        let mut index: HighsInt = rt.scaled.A.mat.start[var as usize];
        while index < rt.scaled.A.mat.start[(var + 1) as usize] {
            if rt.scaled.A.mat.value[index as usize].abs()
                > *maxabscoefperrow
                    .entry(rt.scaled.A.mat.index[index as usize])
                    .or_insert(0.0)
            {
                *maxabscoefperrow
                    .entry(rt.scaled.A.mat.index[index as usize])
                    .or_insert(0.0) = rt.scaled.A.mat.value[index as usize].abs();
            }
            index += 1;
        }
        var += 1;
    }

    let mut var: HighsInt = 0;
    while var < rt.scaled.num_var {
        let mut index: HighsInt = rt.scaled.A.mat.start[var as usize];
        while index < rt.scaled.A.mat.start[(var + 1) as usize] {
            let factor = largestpoweroftwo(
                *maxabscoefperrow
                    .entry(rt.scaled.A.mat.index[index as usize])
                    .or_insert(0.0),
            );
            rt.scaled.A.mat.value[index as usize] /= factor;
            index += 1;
        }
        var += 1;
    }

    let mut row: HighsInt = 0;
    while row < rt.scaled.num_con {
        let factor = largestpoweroftwo(*maxabscoefperrow.entry(row).or_insert(0.0));
        if rt.scaled.con_lo[row as usize] > -f64::INFINITY {
            rt.scaled.con_lo[row as usize] /= factor;
        }
        if rt.scaled.con_up[row as usize] < f64::INFINITY {
            rt.scaled.con_up[row as usize] /= factor;
        }
        row += 1;
    }
}

/// HiGHS: scaling.cpp::scale_cols (lines 60-116)
fn scale_cols(rt: &mut Runtime) {
    if !rt.settings.varscaling {
        return;
    }

    // std::map operator[] -> `*map.entry(key).or_insert(0.0)`.
    let mut maxabscoefpervar: BTreeMap<HighsInt, f64> = BTreeMap::new();
    let mut var: HighsInt = 0;
    while var < rt.scaled.num_var {
        *maxabscoefpervar.entry(var).or_insert(0.0) = 0.0;
        var += 1;
    }

    let mut var: HighsInt = 0;
    while var < rt.scaled.num_var {
        // gather information about variable
        let mut index: HighsInt = rt.scaled.A.mat.start[var as usize];
        while index < rt.scaled.A.mat.start[(var + 1) as usize] {
            if rt.scaled.A.mat.value[index as usize].abs()
                > *maxabscoefpervar.entry(var).or_insert(0.0)
            {
                *maxabscoefpervar.entry(var).or_insert(0.0) =
                    rt.scaled.A.mat.value[index as usize].abs();
            }
            index += 1;
        }

        let mut index: HighsInt = rt.scaled.Q.mat.start[var as usize];
        while index < rt.scaled.Q.mat.start[(var + 1) as usize] {
            if rt.scaled.Q.mat.index[index as usize] == var {
                // PORT NOTE: C fmax
                let cur = *maxabscoefpervar.entry(var).or_insert(0.0);
                *maxabscoefpervar.entry(var).or_insert(0.0) =
                    cur.max(rt.scaled.Q.mat.value[index as usize].abs().sqrt());
            }
            index += 1;
        }
        var += 1;
    }

    let mut factorpervar: BTreeMap<HighsInt, f64> = BTreeMap::new();

    let mut var: HighsInt = 0;
    while var < rt.scaled.num_var {
        // PORT NOTE: C fmin
        let v = 128.0f64.min(largestpoweroftwo(
            *maxabscoefpervar.entry(var).or_insert(0.0),
        ));
        *factorpervar.entry(var).or_insert(0.0) = v;
        var += 1;
    }

    let mut var: HighsInt = 0;
    while var < rt.scaled.num_var {
        // scale variable: constraint matrix
        let factor = *factorpervar.entry(var).or_insert(0.0);
        let mut index: HighsInt = rt.scaled.A.mat.start[var as usize];
        while index < rt.scaled.A.mat.start[(var + 1) as usize] {
            rt.scaled.A.mat.value[index as usize] /= factor;
            index += 1;
        }

        // scale variable: linear objective
        rt.scaled.c.value[var as usize] /= factor;
        var += 1;
    }

    // scale variable: hessian matrix
    let mut var1: HighsInt = 0;
    while var1 < rt.scaled.num_var {
        let factor1 = *factorpervar.entry(var1).or_insert(0.0);
        let mut index: HighsInt = rt.scaled.Q.mat.start[var1 as usize];
        while index < rt.scaled.Q.mat.start[(var1 + 1) as usize] {
            let var2 = rt.scaled.Q.mat.index[index as usize];
            let factor2 = *factorpervar.entry(var2).or_insert(0.0);
            rt.scaled.Q.mat.value[index as usize] /= factor1 * factor2;
            index += 1;
        }
        var1 += 1;
    }
}

/// HiGHS: scaling.cpp::scale (lines 118-123)
pub fn scale(rt: &mut Runtime) {
    rt.scaled = rt.instance.clone();
    scale_rows(rt);
    scale_cols(rt);
    scale_rows(rt);
}

// ---------------------------------------------------------------------------
// Unit tests. Expected values are hand-derived from the C++; if one
// contradicts the C++ semantics, the C++ wins: fix the test and say why
// in the commit.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use crate::qpsolver::*;

    /// min 0.5 (2 x0^2 + 4 x1^2) - 2 x0 - 4 x1, 0 <= x <= 10 and, when
    /// `with_row`, x0 + x1 <= 10. Optimum x = (1, 1), objective -3.
    fn box_qp(with_row: bool) -> Instance {
        let nc = if with_row { 1 } else { 0 };
        let mut inst = Instance::new(2, nc);
        inst.c = QpVector {
            num_nz: 2,
            dim: 2,
            index: vec![0, 1],
            value: vec![-2.0, -4.0],
        };
        inst.Q.mat = MatrixBase {
            num_row: 2,
            num_col: 2,
            start: vec![0, 1, 2],
            index: vec![0, 1],
            value: vec![2.0, 4.0],
        };
        inst.A.mat = if with_row {
            MatrixBase {
                num_row: 1,
                num_col: 2,
                start: vec![0, 1, 2],
                index: vec![0, 0],
                value: vec![1.0, 1.0],
            }
        } else {
            MatrixBase {
                num_row: 0,
                num_col: 2,
                start: vec![0, 0, 0],
                index: vec![],
                value: vec![],
            }
        };
        inst.con_lo = vec![-kHighsInf; nc as usize];
        inst.con_up = vec![10.0; nc as usize];
        inst.var_lo = vec![0.0, 0.0];
        inst.var_up = vec![10.0, 10.0];
        inst
    }

    /// Both variables active at their lower bound (constraint ids
    /// `num_con + i`); the row (id 0) is not in the basis.
    fn at_lower(rt: &mut Runtime) -> Basis {
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kActiveAtLower,
            BasisStatus::kActiveAtLower,
        ];
        Basis::new(rt, vec![1, 2], status, vec![])
    }

    /// Both variables in the nullspace (inactive in basis): Z = I.
    fn all_free(rt: &mut Runtime) -> Basis {
        let status = vec![
            BasisStatus::kInactive,
            BasisStatus::kInactiveInBasis,
            BasisStatus::kInactiveInBasis,
        ];
        Basis::new(rt, vec![], status, vec![1, 2])
    }

    fn close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= 1e-9)
    }

    #[test]
    fn scaling_off_copies_instance() {
        let mut rt = Runtime::new(&box_qp(true), Statistics::default());
        rt.settings.rowscaling = false;
        rt.settings.varscaling = false;
        scale(&mut rt);
        assert_eq!(rt.scaled.c.value, rt.instance.c.value);
        assert_eq!(rt.scaled.A.mat.value, rt.instance.A.mat.value);
    }

    #[test]
    fn largest_power_of_two() {
        // C++ scaling.cpp:16 is `powf(1.0, il)`, which is always 1.0 (the
        // skeleton's expected 4.0 / 0.25 contradicted the C++, so the C++ wins).
        assert_eq!(largestpoweroftwo(5.0), 1.0);
        assert_eq!(largestpoweroftwo(0.3), 1.0);
    }
}
