//! HiGHS source: mip/HighsLpAggregator.h, mip/HighsLpAggregator.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsLpAggregator.h::HighsLpAggregator (lines 26-48)
///
/// Removed C++ members:
/// * C++: `const HighsLpRelaxation & lprelaxation`: removed; methods that use it take
///   `lprelaxation: &HighsLpRelaxation` (and `ms` for `lprelaxation.getRow` /
///   `getMipSolver()`).
#[derive(Clone, Default, Debug)]
pub struct HighsLpAggregator {
    /// C++: `HighsSparseVectorSum vectorsum`
    pub vectorsum: HighsSparseVectorSum,
}

impl HighsLpAggregator {
    /// HiGHS: HighsLpAggregator.cpp::HighsLpAggregator::HighsLpAggregator (lines 13-17)
    /// C++: `HighsLpAggregator(const HighsLpRelaxation & lprelaxation)`
    pub fn new(lprelaxation: &HighsLpRelaxation) -> Self {
        let mut this = HighsLpAggregator::default();
        this.vectorsum
            .set_dimension(lprelaxation.get_lp().num_row_ + lprelaxation.get_lp().num_col_);
        this
    }

    /// HiGHS: HighsLpAggregator.cpp::HighsLpAggregator::addRow (lines 19-28)
    /// C++: `void addRow(HighsInt row, double weight)`
    pub fn add_row(
        &mut self,
        ms: &HighsMipSolver,
        lprelaxation: &HighsLpRelaxation,
        row: HighsInt,
        weight: f64,
    ) {
        let mut len: HighsInt = 0;
        let mut vals: &[f64] = &[];
        let mut inds: &[HighsInt] = &[];
        lprelaxation.get_row(ms, row, &mut len, &mut inds, &mut vals);

        for i in 0..len as usize {
            self.vectorsum.add(inds[i], weight * vals[i]);
        }

        self.vectorsum
            .add(lprelaxation.get_lp().num_col_ + row, -weight);
    }

    /// HiGHS: HighsLpAggregator.cpp::HighsLpAggregator::getCurrentAggregation (lines 30-54)
    /// HiGHS 1.16 (#3507): HighsLpAggregator.cpp::getCurrentAggregation (lines 31-57)
    /// C++: `void getCurrentAggregation(std::vector<HighsInt> & inds, std::vector<double> & vals, bool negate)`
    pub fn get_current_aggregation(
        &mut self,
        ms: &HighsMipSolver,
        lprelaxation: &HighsLpRelaxation,
        inds: &mut Vec<HighsInt>,
        vals: &mut Vec<f64>,
        negate: bool,
    ) {
        let droptol = ms.options_mip_.small_matrix_value;
        // HiGHS 1.16 #3507 (h16_agg_term_zero): HighsLpAggregator.cpp lines 36, 45-46 (1.16.0)
        let h16_agg_term_zero = ms.options_mip_.h16_agg_term_zero;
        let domain: Option<&HighsDomain> = if h16_agg_term_zero {
            Some(ms.get_domain_2())
        } else {
            None
        };
        let num_col = lprelaxation.num_cols();
        self.vectorsum.cleanup(|col: HighsInt, val: f64| {
            // only drop values for columns and not for slack variables of rows as the
            // former are the only ones which might be subject to numerical error. The
            // values for row slack variables are exact copies of the weights used
            // (assuming that separators only add each row once). Also all the
            // separators will only allow rows that do have meaningful coefficient
            // contributions with the used weight
            if let Some(domain) = domain {
                col < num_col
                    && HighsDomain::term_is_zero(
                        val,
                        domain.col_lower_[col as usize],
                        domain.col_upper_[col as usize],
                        droptol,
                    )
            } else {
                col < num_col && val.abs() <= droptol
            }
        });

        inds.clear();
        inds.extend_from_slice(self.vectorsum.get_nonzeros());
        let len = inds.len();
        vals.resize(len, 0.0);

        if negate {
            for i in 0..len {
                vals[i] = -self.vectorsum.get_value(inds[i]);
            }
        } else {
            for i in 0..len {
                vals[i] = self.vectorsum.get_value(inds[i]);
            }
        }
    }

    /// HiGHS: HighsLpAggregator.cpp::HighsLpAggregator::clear (lines 56-56)
    /// C++: `void clear()`
    pub fn clear(&mut self) {
        self.vectorsum.clear();
    }

    /// HiGHS: HighsLpAggregator.h::HighsLpAggregator::isEmpty (lines 47-47)
    /// C++: `bool isEmpty()`
    pub fn is_empty(&mut self) -> bool {
        self.vectorsum.nonzeroinds.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Aggregates two column terms of value 1e-10 (droptol 1e-9): column 0 has bounds
    /// [0, 1e6], column 1 has bounds [0, 1]. Returns the surviving indices and values of
    /// `get_current_aggregation` with `h16_agg_term_zero` set to `on`.
    fn agg_term_zero_run(on: bool) -> (Vec<HighsInt>, Vec<f64>) {
        let mut ms = HighsMipSolver::default();
        {
            let lp = &mut ms.orig_model_;
            lp.num_col_ = 2;
            lp.num_row_ = 0;
            lp.col_cost_ = vec![0.0; 2];
            lp.col_lower_ = vec![0.0, 0.0];
            lp.col_upper_ = vec![1e6, 1.0];
            lp.integrality_ = vec![HighsVarType::kContinuous; 2];
            lp.a_matrix_.num_col_ = 2;
            lp.a_matrix_.start_ = vec![0, 0, 0];
        }
        ms.mipdata_.a_rstart_ = vec![0];
        ms.mipdata_.feastol = 1e-6;
        ms.mipdata_.epsilon = 1e-9;
        ms.options_mip_.small_matrix_value = 1e-9;
        ms.options_mip_.h16_agg_term_zero = on;
        let dom = HighsDomain::new(&ms);
        ms.mipdata_.domains = vec![Slot::new(dom)];

        let mut lprel = HighsLpRelaxation::default();
        lprel.lpsolver.model_lp_.num_col_ = 2;

        let mut aggregator = HighsLpAggregator::default();
        aggregator.vectorsum.set_dimension(2);
        aggregator.vectorsum.add(0, 1e-10);
        aggregator.vectorsum.add(1, 1e-10);

        let (mut inds, mut vals) = (Vec::new(), Vec::new());
        aggregator.get_current_aggregation(&ms, &lprel, &mut inds, &mut vals, false);
        (inds, vals)
    }

    #[test]
    fn h16_agg_term_zero_on_keeps_term_with_large_bound() {
        // 1.16 termIsZero: |1e-10| * max(|0|, |1e6|) = 1e-4 > 1e-9, so column 0 is kept;
        // column 1 has |val| * 1 = 1e-10 <= 1e-9 and is still dropped.
        let (inds, vals) = agg_term_zero_run(true);
        assert_eq!(inds, vec![0]);
        assert_eq!(vals, vec![1e-10]);
    }

    #[test]
    fn h16_agg_term_zero_off_drops_tiny_term() {
        // 1.15.1: every |val| <= droptol is dropped, whatever the bounds.
        let (inds, vals) = agg_term_zero_run(false);
        assert!(inds.is_empty());
        assert!(vals.is_empty());
    }
}
