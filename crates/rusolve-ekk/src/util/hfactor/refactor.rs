//! HiGHS source: util/HFactorRefactor.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Re-factorization from a recorded pivot sequence (`refactor_info_`).

use super::*;
use crate::lp_data::hstruct::RefactorInfo;
use crate::prelude::*;

impl RefactorInfo {
    /// HiGHS: HFactorRefactor.cpp::RefactorInfo::clear (lines 22-28, 7 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    pub fn clear(&mut self) {
        self.r#use = false;
        self.build_synthetic_tick = 0.0;
        self.pivot_var.clear();
        self.pivot_row.clear();
        self.pivot_type.clear();
    }
}

impl HFactor {
    /// HiGHS: HFactorRefactor.cpp::HFactor::rebuild (lines 30-304, 275 lines)
    /// MUST-PORT (determinism, synthetic ticks drive reinversion).
    /// `basic_index` (C++ member pointer `HFactor::basic_index`) is passed explicitly; see `HFactor` docs.
    pub fn rebuild(
        &mut self,
        basic_index: &mut [HighsInt],
        factor_timer_clock_pointer: Option<&HighsTimerClock>,
    ) -> HighsInt {
        let _ = factor_timer_clock_pointer;
        let report_lu = false;
        // Check that the refactorization information should be used
        debug_assert!(self.refactor_info_.r#use);
        // 0. Clear L and U factor
        self.lu_clear();

        self.nwork = 0;
        self.basis_matrix_num_el = 0;
        let mut stage: HighsInt = self.num_row;
        let rank_deficiency: HighsInt;
        let mut has_pivot: Vec<bool> = vec![false; self.num_row as usize];
        let report_unit = false;
        let report_singletons = false;
        let report_markowitz = false;
        let report_anything = report_unit || report_singletons || report_markowitz;
        if report_anything {
            println!("\nRefactor");
        }
        // Take build_synthetic_tick from the refactor info so that this
        // refactorization doesn't look unrealistically cheap.
        self.build_synthetic_tick = self.refactor_info_.build_synthetic_tick;
        // Check that the refactorization info has been set up
        debug_assert!(self.refactor_info_.pivot_row.len() as HighsInt >= self.num_row);
        debug_assert!(self.refactor_info_.pivot_var.len() as HighsInt >= self.num_row);
        debug_assert!(self.refactor_info_.pivot_type.len() as HighsInt >= self.num_row);
        for i_k in 0..self.num_row {
            let i_row: HighsInt = self.refactor_info_.pivot_row[i_k as usize];
            let i_var: HighsInt = self.refactor_info_.pivot_var[i_k as usize];
            let pivot_type: i8 = self.refactor_info_.pivot_type[i_k as usize];
            debug_assert!(!has_pivot[i_row as usize]);

            if pivot_type == kPivotLogical || pivot_type == kPivotUnit {
                if pivot_type == kPivotLogical {
                    // 1.1 Logical column
                    if report_unit {
                        println!("Stage {}: Logical", i_k);
                    }
                    debug_assert!(i_var >= self.num_col);
                    self.basis_matrix_num_el += 1;
                } else if pivot_type == kPivotUnit {
                    // 1.2 (Structural) unit column
                    if report_unit {
                        println!("Stage {}: Unit", i_k);
                    }
                    debug_assert!(i_var < self.num_col);
                    let start: HighsInt = self.a_start[i_var as usize];
                    let count: HighsInt = self.a_start[(i_var + 1) as usize] - start;
                    debug_assert!(self.a_index[start as usize] == i_row);
                    debug_assert!(count == 1 && self.a_value[start as usize] == 1.0);
                    let _ = count;
                    self.basis_matrix_num_el += 1;
                }
                // 1.3 Record unit column
                self.l_start.push(self.l_index.len() as HighsInt);
                self.u_pivot_index.push(i_row);
                self.u_pivot_value.push(1.0);
                self.u_start.push(self.u_index.len() as HighsInt);
            } else if pivot_type == kPivotRowSingleton || pivot_type == kPivotColSingleton {
                // Row or column singleton
                debug_assert!(i_var < self.num_col);
                let start: HighsInt = self.a_start[i_var as usize];
                let end: HighsInt = self.a_start[(i_var + 1) as usize];
                // Find where the pivot is
                let mut pivot_k: HighsInt = -1;
                for k in start..end {
                    if self.a_index[k as usize] == i_row {
                        pivot_k = k;
                        break;
                    }
                }
                debug_assert!(pivot_k >= 0);
                // Check that the pivot isn't too small. Shouldn't happen since
                // this is refactorization
                let abs_pivot: f64 = self.a_value[pivot_k as usize].abs();
                debug_assert!(abs_pivot >= self.pivot_tolerance);
                if abs_pivot < self.pivot_tolerance {
                    rank_deficiency = self.nwork + 1;
                    return rank_deficiency;
                }
                if pivot_type == kPivotRowSingleton {
                    // 2.2 Deal with row singleton
                    let pivot_multiplier: f64 = 1.0 / self.a_value[pivot_k as usize];
                    if report_singletons {
                        println!(
                            "Stage {}: Row singleton ({:4}, {})",
                            i_k, pivot_k, pivot_multiplier
                        );
                    }
                    for section in 0..2 {
                        let p0: HighsInt = if section == 0 { start } else { pivot_k + 1 };
                        let p1: HighsInt = if section == 0 { pivot_k } else { end };
                        for k in p0..p1 {
                            let local_i_row: HighsInt = self.a_index[k as usize];
                            if !has_pivot[local_i_row as usize] {
                                if report_singletons {
                                    println!(
                                        "Row singleton: L En ({:4}, {:11.4})",
                                        local_i_row,
                                        self.a_value[k as usize] * pivot_multiplier
                                    );
                                }
                                self.l_index.push(local_i_row);
                                self.l_value
                                    .push(self.a_value[k as usize] * pivot_multiplier);
                            } else {
                                if report_singletons {
                                    println!(
                                        "Row singleton: U En ({:4}, {:11.4})",
                                        local_i_row, self.a_value[k as usize]
                                    );
                                }
                                self.u_index.push(local_i_row);
                                self.u_value.push(self.a_value[k as usize]);
                            }
                        }
                    }
                    self.l_start.push(self.l_index.len() as HighsInt);
                    if report_singletons {
                        println!(
                            "Row singleton: U Pv ({:4}, {:11.4})",
                            i_row, self.a_value[pivot_k as usize]
                        );
                    }
                    self.u_pivot_index.push(i_row);
                    self.u_pivot_value.push(self.a_value[pivot_k as usize]);
                    self.u_start.push(self.u_index.len() as HighsInt);
                } else {
                    // 2.3 Deal with column singleton
                    if report_singletons {
                        println!("Stage {}: Col singleton", i_k);
                    }
                    for k in start..pivot_k {
                        if report_singletons {
                            println!(
                                "Col singleton: U En ({:4}, {:11.4})",
                                self.a_index[k as usize], self.a_value[k as usize]
                            );
                        }
                        self.u_index.push(self.a_index[k as usize]);
                        self.u_value.push(self.a_value[k as usize]);
                    }
                    for k in (pivot_k + 1)..end {
                        if report_singletons {
                            println!(
                                "Col singleton: U En ({:4}, {:11.4})",
                                self.a_index[k as usize], self.a_value[k as usize]
                            );
                        }
                        self.u_index.push(self.a_index[k as usize]);
                        self.u_value.push(self.a_value[k as usize]);
                    }
                    self.l_start.push(self.l_index.len() as HighsInt);
                    if report_singletons {
                        println!(
                            "Col singleton: U Pv ({:4}, {:11.4})",
                            i_row, self.a_value[pivot_k as usize]
                        );
                    }
                    self.u_pivot_index.push(i_row);
                    self.u_pivot_value.push(self.a_value[pivot_k as usize]);
                    self.u_start.push(self.u_index.len() as HighsInt);
                }
            } else {
                debug_assert!(pivot_type == kPivotMarkowitz);
                stage = i_k;
                break;
            }
            basic_index[i_row as usize] = i_var;
            has_pivot[i_row as usize] = true;
        }
        if report_lu {
            println!("\nAfter units and singletons");
            self.report_lu(kReportLuBoth, false);
        }
        if stage < self.num_row {
            // Handle the remaining Markowitz pivots
            //
            // First of all complete the L factor with identity columns so
            // that FtranL counts the RHS entries in rows that don't yet have
            // pivots by running to completion. In the hyper-sparse code,
            // these will HOPEFULLY be skipped
            //
            // There are already l_start entries for the first stage rows, but
            // l_pivot_index is not assigned, as u_pivot_index gets copied into it
            self.l_start.resize((self.num_row + 1) as usize, 0);
            for i_k in stage..self.num_row {
                self.l_start[(i_k + 1) as usize] = self.l_start[i_k as usize];
            }
            self.l_pivot_index.resize(self.num_row as usize, 0);
            for i_k in 0..self.num_row {
                self.l_pivot_index[i_k as usize] = self.refactor_info_.pivot_row[i_k as usize];
            }
            // To do hyper-sparse FtranL operations, have to set up l_pivot_lookup.
            self.l_pivot_lookup.resize(self.num_row as usize, 0);
            for i_row in 0..self.num_row {
                if i_row < stage
                    && self.l_pivot_lookup[self.l_pivot_index[i_row as usize] as usize] != i_row
                {
                    // printf("Strange: Thought that l_pivot_lookup[l_pivot_index[iRow]] == iRow\n");
                }
                self.l_pivot_lookup[self.l_pivot_index[i_row as usize] as usize] = i_row;
            }
            // Need to know whether to consider matrix entries for FtranL
            // operation. Initially these correspond to all the rows without
            // pivots
            let not_in_bump: Vec<bool> = has_pivot.clone();
            // Monitor density of FtranL result to possibly switch from exploiting
            // hyper-sparsity
            let mut expected_density: f64 = 0.0;
            // Initialise a HVector in which the L and U entries of the
            // pivotal column will be formed
            let mut column = HVector::default();
            column.setup(self.num_row);
            for i_k in stage..self.num_row {
                let i_row: HighsInt = self.refactor_info_.pivot_row[i_k as usize];
                let i_var: HighsInt = self.refactor_info_.pivot_var[i_k as usize];
                let pivot_type: i8 = self.refactor_info_.pivot_type[i_k as usize];
                debug_assert!(!has_pivot[i_row as usize]);
                debug_assert!(pivot_type == kPivotMarkowitz);
                let _ = pivot_type;
                // Set up the column for the FtranL. It contains the matrix
                // entries in rows without pivots, and the remaining entries
                // start forming the U column
                column.clear();
                let mut start: HighsInt = self.a_start[i_var as usize];
                let mut end: HighsInt = self.a_start[(i_var + 1) as usize];
                for i_el in start..end {
                    let local_i_row: HighsInt = self.a_index[i_el as usize];
                    if not_in_bump[local_i_row as usize] {
                        self.u_index.push(local_i_row);
                        self.u_value.push(self.a_value[i_el as usize]);
                    } else {
                        column.index[column.count as usize] = local_i_row;
                        column.count += 1;
                        column.array[local_i_row as usize] = self.a_value[i_el as usize];
                    }
                }
                // Perform FtranL, but don't time it!
                self.ftran_l(&mut column, expected_density, None);
                // Update the running average density
                let local_density: f64 = (1.0 * column.count as f64) / self.num_row as f64;
                expected_density = kRunningAverageMultiplier * local_density
                    + (1.0 - kRunningAverageMultiplier) * expected_density;
                // Strip out small values
                column.tight();
                // Now form the column of L
                //
                // Find the pivot
                let mut pivot_k: HighsInt = -1;
                start = 0;
                end = column.count;
                for k in start..end {
                    if column.index[k as usize] == i_row {
                        pivot_k = k;
                        break;
                    }
                }
                debug_assert!(pivot_k >= 0);
                // Check that the pivot isn't too small. Shouldn't happen since
                // this is refactorization
                let abs_pivot: f64 = column.array[i_row as usize].abs();
                debug_assert!(abs_pivot >= self.pivot_tolerance);
                if abs_pivot < self.pivot_tolerance {
                    rank_deficiency = self.num_row - i_k;
                    return rank_deficiency;
                }
                let pivot_multiplier: f64 = 1.0 / column.array[i_row as usize];
                for section in 0..2 {
                    let p0: HighsInt = if section == 0 { start } else { pivot_k + 1 };
                    let p1: HighsInt = if section == 0 { pivot_k } else { end };
                    for k in p0..p1 {
                        let local_i_row: HighsInt = column.index[k as usize];
                        if !has_pivot[local_i_row as usize] {
                            self.l_index.push(local_i_row);
                            self.l_value
                                .push(column.array[local_i_row as usize] * pivot_multiplier);
                        } else {
                            self.u_index.push(local_i_row);
                            self.u_value.push(column.array[local_i_row as usize]);
                        }
                    }
                }
                self.l_start[(i_k + 1) as usize] = self.l_index.len() as HighsInt;
                self.u_pivot_index.push(i_row);
                self.u_pivot_value.push(column.array[i_row as usize]);
                self.u_start.push(self.u_index.len() as HighsInt);
                basic_index[i_row as usize] = i_var;
                has_pivot[i_row as usize] = true;
                if report_lu {
                    println!("\nAfter Markowitz {}", i_k - stage);
                    self.report_lu(kReportLuBoth, false);
                }
            }
        }
        if report_lu {
            println!("\nRefactored INVERT");
            self.report_lu(kReportLuBoth, false);
        }
        self.build_finish(basic_index);
        0
    }
}
