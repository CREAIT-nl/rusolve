//! HiGHS source: simplex/HEkk.cpp (lines 4054-4409)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Infeasibility proofs, value scaling, iterative refinement, and the HighsSimplexStats /
//! HighsRayRecord members defined at the end of HEkk.cpp.

use crate::prelude::*;

impl HEkk {
    /// HiGHS: HEkk.cpp::HEkk::proofOfPrimalInfeasibility (lines 4054-4065, 12 lines)
    pub fn proof_of_primal_infeasibility(&mut self) -> bool {
        // To be called from outside HEkk when row_ep is not known
        debug_assert!(self.dual_ray_record_.index >= 0);
        let move_out = self.dual_ray_record_.sign;
        let row_out = self.dual_ray_record_.index;
        // Compute the basis inverse row
        let mut row_ep = HVector::default();
        row_ep.setup(self.lp_.num_row_);
        self.unit_btran(row_out, &mut row_ep);
        self.proof_of_primal_infeasibility_ep(&mut row_ep, move_out, row_out)
    }

    /// HiGHS: HEkk.cpp::HEkk::proofOfPrimalInfeasibility (lines 4067-4260, 194 lines)
    pub fn proof_of_primal_infeasibility_ep(
        &mut self,
        row_ep: &mut HVector,
        move_out: HighsInt,
        row_out: HighsInt,
    ) -> bool {
        // To be called from inside HEkkDual
        let mut debug_product_report: HighsInt = kDebugReportOff;
        let debug_proof_report_on = true;
        let mut debug_rows_report = false;
        let mut debug_proof_report = false;
        if self.debug_iteration_report_ {
            if debug_proof_report_on {
                debug_product_report = kDebugReportOff; // kDebugReportAll; //
            }
            debug_rows_report = debug_proof_report_on;
            debug_proof_report = debug_proof_report_on;
        }

        let use_row_wise_matrix = self.status_.has_ar_matrix;
        let use_iterative_refinement = false; // debug_iteration_report_;//
        if use_iterative_refinement {
            self.simplex_nla_.report_array("Row e_p.0", row_ep, true);
            self.unit_btran_iterative_refinement(row_out, row_ep);
            self.simplex_nla_.report_array("Row e_p.1", row_ep, true);
        }

        // Refine row_ep by removing relatively small values
        let row_ep_scale: f64 = 0.0;
        // if (use_refinement) refineArray(row_ep, row_ep_scale,
        // refinement_tolerance);
        // Determine the maximum absolute value in row_ep
        let mut proof_lower = HighsCDouble::from(0.0);
        let max_num_basic_proof_report: HighsInt = 25;
        let max_num_zeroed_report: HighsInt = 25;
        let mut num_zeroed_for_small_report: HighsInt = 0;
        let mut max_zeroed_for_small_value: f64 = 0.0;
        let mut num_zeroed_for_lb_report: HighsInt = 0;
        let mut max_zeroed_for_lb_value: f64 = 0.0;
        let mut num_zeroed_for_ub_report: HighsInt = 0;
        let mut max_zeroed_for_ub_value: f64 = 0.0;
        for i_x in 0..row_ep.count {
            let i_row = row_ep.index[i_x as usize];
            // Give row_ep the sign of the leaving row - as is done in
            // getDualRayInterface.
            let row_ep_value = row_ep.array[i_row as usize];
            debug_assert!(row_ep_value != 0.0);
            if (row_ep_value * self.get_max_abs_row_value(i_row)).abs()
                <= self.options_.small_matrix_value
            {
                if debug_proof_report {
                    let abs_row_ep_value = row_ep_value.abs();
                    if num_zeroed_for_small_report < max_num_zeroed_report
                        && max_zeroed_for_small_value < abs_row_ep_value
                    {
                        println!(
                            "Zeroed row_ep.array[{:6}] = {:11.4e} due to being small in contribution",
                            i_row, row_ep_value
                        );
                        num_zeroed_for_small_report += 1;
                        max_zeroed_for_small_value = abs_row_ep_value;
                    }
                }
                row_ep.array[i_row as usize] = 0.0;
                continue;
            }

            row_ep.array[i_row as usize] *= move_out as f64;

            // make sure infinite sides are not used
            let row_bound: f64;
            if row_ep.array[i_row as usize] > 0.0 {
                row_bound = self.lp_.row_lower_[i_row as usize];
                if highs_is_infinity(-row_bound) {
                    // row lower bound is infinite
                    if debug_proof_report {
                        let abs_row_ep_value = row_ep_value.abs();
                        if num_zeroed_for_lb_report < max_num_zeroed_report
                            && max_zeroed_for_lb_value < abs_row_ep_value
                        {
                            println!(
                                "Zeroed row_ep.array[{:6}] = {:11.4e} due to infinite lower bound",
                                i_row, row_ep_value
                            );
                            num_zeroed_for_lb_report += 1;
                            max_zeroed_for_lb_value = abs_row_ep_value;
                        }
                    }
                    row_ep.array[i_row as usize] = 0.0;
                    continue;
                }
            } else {
                row_bound = self.lp_.row_upper_[i_row as usize];
                if highs_is_infinity(row_bound) {
                    // row upper bound is infinite
                    if debug_proof_report {
                        let abs_row_ep_value = row_ep_value.abs();
                        if num_zeroed_for_ub_report < max_num_zeroed_report
                            && max_zeroed_for_ub_value < abs_row_ep_value
                        {
                            println!(
                                "Zeroed row_ep.array[{:6}] = {:11.4e} due to infinite upper bound",
                                i_row, row_ep_value
                            );
                            num_zeroed_for_ub_report += 1;
                            max_zeroed_for_ub_value = abs_row_ep_value;
                        }
                    }
                    row_ep.array[i_row as usize] = 0.0;
                    continue;
                }
            }
            // add up lower bound of proof constraint
            proof_lower += row_ep.array[i_row as usize] * row_bound;
        }
        // Form the proof constraint coefficients
        self.proof_value_.clear();
        self.proof_index_.clear();
        if use_row_wise_matrix {
            self.ar_matrix_.product_transpose_quad_hvector(
                &mut self.proof_value_,
                &mut self.proof_index_,
                row_ep,
                debug_product_report,
            );
        } else {
            self.lp_.a_matrix_.product_transpose_quad_hvector(
                &mut self.proof_value_,
                &mut self.proof_index_,
                row_ep,
                debug_product_report,
            );
        }

        let proof_num_nz = self.proof_index_.len() as HighsInt;
        if debug_rows_report {
            self.simplex_nla_.report_array("Row e_p", row_ep, true);
            self.simplex_nla_.report_vector(
                "Proof",
                proof_num_nz,
                &self.proof_value_,
                &self.proof_index_,
                true,
            );
        }
        if debug_proof_report {
            println!(
                "HEkk::proofOfPrimalInfeasibility row_ep.count = {}; proof_num_nz = {}; row_ep_scale = {}",
                row_ep.count, proof_num_nz, row_ep_scale
            );
            let mut num_basic_proof_report: HighsInt = 0;
            let mut max_basic_proof_value: f64 = 0.0;
            for i in 0..proof_num_nz {
                let i_col = self.proof_index_[i as usize];
                let value = self.proof_value_[i as usize];
                let abs_value = value.abs();
                if self.basis_.nonbasic_flag_[i_col as usize] == 0
                    && max_basic_proof_value < abs_value
                    && num_basic_proof_report < max_num_basic_proof_report
                {
                    println!(
                        "Proof entry {:6} (Column {:6}) is basic with value {:11.4e}",
                        i, i_col, value
                    );
                    max_basic_proof_value = abs_value;
                    num_basic_proof_report += 1;
                }
            }
        }
        let mut implied_upper = HighsCDouble::from(0.0);
        let mut sum_inf = HighsCDouble::from(0.0);
        for i in 0..proof_num_nz {
            let i_col = self.proof_index_[i as usize];
            let value = self.proof_value_[i as usize];
            if value > 0.0 {
                if highs_is_infinity(self.lp_.col_upper_[i_col as usize]) {
                    sum_inf += value;
                    if sum_inf > self.options_.small_matrix_value {
                        break;
                    }
                    continue;
                    // Commented out unreachable code
                    //        if (value <= options_->small_matrix_value) continue;
                }
                implied_upper += value * self.lp_.col_upper_[i_col as usize];
            } else {
                if highs_is_infinity(-self.lp_.col_lower_[i_col as usize]) {
                    sum_inf += -value;
                    if sum_inf > self.options_.small_matrix_value {
                        break;
                    }
                    continue;
                }
                implied_upper += value * self.lp_.col_lower_[i_col as usize];
            }
        }
        let infinite_implied_upper = sum_inf > self.options_.small_matrix_value;
        let gap = f64::from(proof_lower - implied_upper);
        let gap_ok = gap > self.options_.primal_feasibility_tolerance;
        let proof_of_primal_infeasibility = !infinite_implied_upper && gap_ok;

        let local_report = false;
        if !proof_of_primal_infeasibility && local_report {
            println!(
                "HEkk::proofOfPrimalInfeasibility: row {:6}; gap = {:11.4e} ({}); sumInf = {:11.4e} ({}) so proof is {}",
                row_out,
                gap,
                highs_bool_to_string(gap_ok, 0),
                f64::from(sum_inf),
                highs_bool_to_string(infinite_implied_upper, 0),
                highs_bool_to_string(proof_of_primal_infeasibility, 0)
            );
        }
        if debug_proof_report {
            print!(
                "HEkk::proofOfPrimalInfeasibility has {}finite implied upper bound",
                if infinite_implied_upper { "in" } else { "" }
            );
            if !infinite_implied_upper {
                print!(" and gap = {}", gap);
            }
            println!(
                " so proof is {}",
                if proof_of_primal_infeasibility {
                    "true"
                } else {
                    "false"
                }
            );
        }
        proof_of_primal_infeasibility
    }

    /// HiGHS: HEkk.cpp::HEkk::getValueScale (lines 4262-4269, 8 lines)
    pub fn get_value_scale(&self, count: HighsInt, value: &[f64]) -> f64 {
        if count <= 0 {
            return 1.0;
        }
        let mut max_abs_value: f64 = 0.0;
        for i_x in 0..count {
            max_abs_value = cpp_max(value[i_x as usize].abs(), max_abs_value);
        }
        nearest_power_of_two_scale(max_abs_value)
    }

    /// HiGHS: HEkk.cpp::HEkk::getMaxAbsRowValue (lines 4271-4279, 9 lines)
    pub fn get_max_abs_row_value(&mut self, row: HighsInt) -> f64 {
        if !self.status_.has_ar_matrix {
            self.initialise_partitioned_rowwise_matrix();
        }

        let mut val: f64 = -1.0;
        for i in self.ar_matrix_.start_[row as usize]..self.ar_matrix_.start_[(row + 1) as usize] {
            val = cpp_max(val, self.ar_matrix_.value_[i as usize].abs());
        }

        val
    }

    /// HiGHS: HEkk.cpp::HEkk::unitBtranIterativeRefinement (lines 4281-4328, 48 lines)
    pub fn unit_btran_iterative_refinement(&mut self, row_out: HighsInt, row_ep: &mut HVector) {
        // Perform an iteration of refinement
        let mut residual = HVector::default();
        let mut residual_norm: f64 = 0.0;
        let mut correction_norm: f64;
        let expected_density: f64 = 1.0;
        residual.setup(self.lp_.num_row_);
        self.unit_btran_residual(row_out, row_ep, &mut residual, &mut residual_norm);
        let debug_iterative_refinement_report_on = false;
        let mut debug_iterative_refinement_report = false;
        if self.debug_iteration_report_ {
            debug_iterative_refinement_report = debug_iterative_refinement_report_on;
        }
        if debug_iterative_refinement_report {
            println!(
                "HEkk::unitBtranIterativeRefinement: Residual   has {:6} / {:6} nonzeros and norm of {}",
                residual.count, self.lp_.num_row_, residual_norm
            );
        }
        if residual_norm == 0.0 {
            return;
        }
        // Normalise using nearest power of 2 to ||correction_rhs|| so kHighsTiny
        // isn't used adversely
        let residual_scale = nearest_power_of_two_scale(residual_norm);
        for i_el in 0..residual.count {
            let idx = residual.index[i_el as usize] as usize;
            residual.array[idx] *= residual_scale;
        }
        self.btran(&mut residual, expected_density);
        row_ep.count = 0;
        correction_norm = 0.0;
        // Adding two (possibly sparse) vectors, so have to loop over all rows
        for i_row in 0..self.lp_.num_row_ {
            if residual.array[i_row as usize] != 0.0 {
                let correction_value = residual.array[i_row as usize] / residual_scale;
                correction_norm = cpp_max(correction_value.abs(), correction_norm);
                row_ep.array[i_row as usize] -= correction_value;
            }
            if row_ep.array[i_row as usize].abs() < kHighsTiny {
                row_ep.array[i_row as usize] = 0.0;
            } else {
                row_ep.index[row_ep.count as usize] = i_row;
                row_ep.count += 1;
            }
        }
        if debug_iterative_refinement_report {
            println!(
                "HEkk::unitBtranIterativeRefinement: Correction has {:6} / {:6} nonzeros and norm of {}",
                residual.count, self.lp_.num_row_, correction_norm
            );
        }
    }

    /// HiGHS: HEkk.cpp::HEkk::unitBtranResidual (lines 4330-4360, 31 lines)
    pub fn unit_btran_residual(
        &mut self,
        row_out: HighsInt,
        row_ep: &HVector,
        residual: &mut HVector,
        residual_norm: &mut f64,
    ) {
        let lp = &self.lp_;
        let mut quad_residual: Vec<HighsCDouble> =
            vec![HighsCDouble::from(0.0); lp.num_row_ as usize];
        quad_residual[row_out as usize] = HighsCDouble::from(-1.0);
        for i_row in 0..lp.num_row_ {
            let i_var = self.basis_.basic_index_[i_row as usize];
            let mut value = quad_residual[i_row as usize];
            if i_var < lp.num_col_ {
                for i_el in
                    lp.a_matrix_.start_[i_var as usize]..lp.a_matrix_.start_[(i_var + 1) as usize]
                {
                    value += lp.a_matrix_.value_[i_el as usize]
                        * row_ep.array[lp.a_matrix_.index_[i_el as usize] as usize];
                }
            } else {
                value += row_ep.array[(i_var - lp.num_col_) as usize];
            }
            quad_residual[i_row as usize] = value;
        }
        residual.clear();
        residual.pack_flag = false;
        *residual_norm = 0.0;
        for i_row in 0..lp.num_row_ {
            let value = f64::from(quad_residual[i_row as usize]);
            if value != 0.0 {
                residual.array[i_row as usize] = value;
                residual.index[residual.count as usize] = i_row;
                residual.count += 1;
            }
            *residual_norm = cpp_max(residual.array[i_row as usize].abs(), *residual_norm);
        }
    }
}

impl HighsSimplexStats {
    /// HiGHS: HEkk.cpp::HighsSimplexStats::report (lines 4362-4377, 16 lines)
    /// Reporting: prints the stats; keep (used for traces).
    pub fn report(&self, file: Option<HighsLogStream>, message: &str) {
        let Some(_file) = file else {
            return;
        };
        println!("\nSimplex stats: {}", message);
        println!("   valid                      = {}", self.valid as i32);
        println!("   iteration_count            = {}", self.iteration_count);
        println!("   num_invert                 = {}", self.num_invert);
        println!(
            "   last_invert_num_el         = {}",
            self.last_invert_num_el
        );
        println!(
            "   last_factored_basis_num_el = {}",
            self.last_factored_basis_num_el
        );
        println!("   col_aq_density             = {}", self.col_aq_density);
        println!("   row_ep_density             = {}", self.row_ep_density);
        println!("   row_ap_density             = {}", self.row_ap_density);
        println!("   row_DSE_density            = {}", self.row_dse_density);
    }

    /// HiGHS: HEkk.cpp::HighsSimplexStats::initialise (lines 4379-4389, 11 lines)
    pub fn initialise(&mut self, iteration_count_: HighsInt) {
        self.valid = false;
        self.iteration_count = -iteration_count_;
        self.num_invert = 0;
        self.last_invert_num_el = 0;
        self.last_factored_basis_num_el = 0;
        self.col_aq_density = 0.0;
        self.row_ep_density = 0.0;
        self.row_ap_density = 0.0;
        self.row_dse_density = 0.0;
    }
}

impl HighsRayRecord {
    /// HiGHS: HEkk.cpp::HighsRayRecord::getRayRecord (lines 4391-4397, 7 lines)
    pub fn get_ray_record(&self) -> HighsRayRecord {
        HighsRayRecord {
            index: self.index,
            sign: self.sign,
            value: self.value.clone(),
        }
    }

    /// HiGHS: HEkk.cpp::HighsRayRecord::setRayRecord (lines 4399-4403, 5 lines)
    pub fn set_ray_record(&mut self, from_record: &HighsRayRecord) {
        self.index = from_record.index;
        self.sign = from_record.sign;
        self.value = from_record.value.clone();
    }

    /// HiGHS: HEkk.cpp::HighsRayRecord::clear (lines 4405-4409, 5 lines)
    pub fn clear(&mut self) {
        self.index = kNoRayIndex;
        self.sign = kNoRaySign;
        self.value.clear();
    }
}
