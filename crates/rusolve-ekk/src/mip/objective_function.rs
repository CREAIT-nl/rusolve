//! HiGHS source: mip/HighsObjectiveFunction.h, mip/HighsObjectiveFunction.cpp
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Part of the HiGHS MIP port.

use crate::prelude::*;

/// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction (lines 22-69)
///
/// C++ `const HighsLp * model` is kept as the selector `model` of the model it points at
/// (`mipsolver.model_` at construction); methods that read it take `model: &HighsLp`, and
/// callers pass the model `model` selects.
#[derive(Clone, Default, Debug)]
pub struct HighsObjectiveFunction {
    /// C++: `const HighsLp * model`
    pub model: HighsMipModel,
    /// C++: `double objIntScale`
    pub obj_int_scale: f64,
    /// C++: `HighsInt numIntegral`
    pub num_integral: HighsInt,
    /// C++: `HighsInt numBinary`
    pub num_binary: HighsInt,
    /// C++: `std::vector<HighsInt> objectiveNonzeros`
    pub objective_nonzeros: Vec<HighsInt>,
    /// C++: `std::vector<double> objectiveVals`
    pub objective_vals: Vec<f64>,
    /// C++: `std::vector<HighsInt> cliquePartitionStart`
    pub clique_partition_start: Vec<HighsInt>,
    /// C++: `std::vector<HighsInt> colToPartition`
    pub col_to_partition: Vec<HighsInt>,
}

impl HighsObjectiveFunction {
    /// HiGHS: HighsObjectiveFunction.cpp::HighsObjectiveFunction::HighsObjectiveFunction (lines 20-62)
    /// C++: `HighsObjectiveFunction(const HighsMipSolver & mipsolver)`
    pub fn new(mipsolver: &HighsMipSolver) -> Self {
        let model = crate::mip_model!(mipsolver);

        let mut objective_nonzeros: Vec<HighsInt> = Vec::with_capacity(model.num_col_ as usize);
        for i in 0..model.num_col_ {
            if model.col_cost_[i as usize] != 0.0 {
                objective_nonzeros.push(i);
            }
        }

        let col_to_partition = vec![-1; model.num_col_ as usize];
        let mut this = Self {
            model: mipsolver.model_,
            obj_int_scale: 0.0,
            num_integral: 0,
            num_binary: 0,
            objective_nonzeros,
            objective_vals: Vec::new(),
            clique_partition_start: vec![0],
            col_to_partition,
        };

        if this.objective_nonzeros.is_empty() {
            this.num_binary = 0;
            this.num_integral = 0;
            this.obj_int_scale = 1.0;
            return this;
        }

        this.num_integral = cpp_partition(&mut this.objective_nonzeros[..], |i: &HighsInt| {
            !mipsolver.is_col_continuous(*i)
        }) as HighsInt;

        if this.num_integral == 0 {
            this.num_binary = 0;
        } else {
            let num_integral = this.num_integral as usize;
            this.num_binary = cpp_partition(
                &mut this.objective_nonzeros[..num_integral],
                |i: &HighsInt| {
                    model.col_lower_[*i as usize] == 0.0 && model.col_upper_[*i as usize] == 1.0
                },
            ) as HighsInt;
        }

        for &i in this.objective_nonzeros.iter() {
            this.objective_vals.push(model.col_cost_[i as usize]);
        }

        this.obj_int_scale = 0.0;
        this
    }

    /// HiGHS: HighsObjectiveFunction.cpp::HighsObjectiveFunction::setupCliquePartition (lines 64-113)
    /// C++: `void setupCliquePartition(const HighsDomain & globaldom, HighsCliqueTable & cliqueTable)`
    pub fn setup_clique_partition(
        &mut self,
        model: &HighsLp,
        globaldom: &HighsDomain,
        clique_table: &mut HighsCliqueTable,
    ) {
        // skip if not more than 1 binary column is present
        if self.num_binary <= 1 {
            return;
        }

        let mut clqvars: Vec<HighsCliqueTable_CliqueVar> = Vec::new();

        let binary_end = self.num_binary as usize;
        for it in 0..binary_end {
            let col = self.objective_nonzeros[it];
            clqvars.push(HighsCliqueTable_CliqueVar::new(
                col,
                (model.col_cost_[col as usize] < 0.0) as HighsInt,
            ));
        }

        clique_table.clique_partition_2(
            &model.col_cost_,
            &mut clqvars,
            &mut self.clique_partition_start,
        );
        let num_partitions = self.clique_partition_start.len() as HighsInt - 1;
        if num_partitions == self.num_binary {
            self.clique_partition_start.resize(1, 0);
        } else {
            let mut p: HighsInt = 0;
            let mut k: HighsInt = 0;

            for i in 0..num_partitions {
                let iu = i as usize;
                if self.clique_partition_start[iu + 1] - self.clique_partition_start[iu] == 1 {
                    continue;
                }

                self.clique_partition_start[p as usize] = k;

                let mut j = self.clique_partition_start[iu];
                while j < self.clique_partition_start[iu + 1] {
                    self.col_to_partition[clqvars[j as usize].col() as usize] = k;
                    k += 1;
                    j += 1;
                }

                p += 1;
            }

            self.clique_partition_start[p as usize] = k;
            self.clique_partition_start.resize((p + 1) as usize, 0);

            let col_to_partition = &self.col_to_partition;
            let num_binary = self.num_binary as usize;
            crate::util::highs_pdqsort::pdqsort_by(
                &mut self.objective_nonzeros[..num_binary],
                |i: &HighsInt, j: &HighsInt| {
                    (
                        col_to_partition[*i as usize] as HighsUInt,
                        crate::util::highs_hash::hash_i32(*i),
                    ) < (
                        col_to_partition[*j as usize] as HighsUInt,
                        crate::util::highs_hash::hash_i32(*j),
                    )
                },
            );

            for i in 0..num_binary {
                self.objective_vals[i] = model.col_cost_[self.objective_nonzeros[i] as usize];
            }
        }
    }

    /// HiGHS: HighsObjectiveFunction.cpp::HighsObjectiveFunction::checkIntegrality (lines 115-124)
    /// C++: `void checkIntegrality(double epsilon)`
    pub fn check_integrality(&mut self, epsilon: f64) {
        if self.num_integral == self.objective_nonzeros.len() as HighsInt {
            if self.num_integral != 0 {
                self.obj_int_scale =
                    HighsIntegers::integral_scale_vec(&self.objective_vals, epsilon, epsilon);
                if self.obj_int_scale * kHighsTiny > epsilon {
                    self.obj_int_scale = 0.0;
                }
            } else {
                self.obj_int_scale = 1.0;
            }
        }
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::getObjectiveNonzeros (lines 42-44)
    /// C++: `const std::vector<HighsInt> & getObjectiveNonzeros() const`
    pub fn get_objective_nonzeros(&self) -> &[HighsInt] {
        &self.objective_nonzeros
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::getObjectiveValuesPacked (lines 46-48)
    /// C++: `const std::vector<double> & getObjectiveValuesPacked() const`
    pub fn get_objective_values_packed(&self) -> &[f64] {
        &self.objective_vals
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::getNumBinariesInObjective (lines 50-50)
    /// C++: `HighsInt getNumBinariesInObjective() const`
    pub fn get_num_binaries_in_objective(&self) -> HighsInt {
        self.num_binary
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::getCliquePartitionStarts (lines 52-54)
    /// C++: `const std::vector<HighsInt> & getCliquePartitionStarts() const`
    pub fn get_clique_partition_starts(&self) -> &[HighsInt] {
        &self.clique_partition_start
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::getNumCliquePartitions (lines 56-58)
    /// C++: `HighsInt getNumCliquePartitions() const`
    pub fn get_num_clique_partitions(&self) -> HighsInt {
        self.clique_partition_start.len() as HighsInt - 1
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::getColCliquePartition (lines 60-62)
    /// C++: `HighsInt getColCliquePartition(HighsInt col) const`
    pub fn get_col_clique_partition(&self, col: HighsInt) -> HighsInt {
        self.col_to_partition[col as usize]
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::integralScale (lines 64-64)
    /// C++: `double integralScale() const`
    pub fn integral_scale(&self) -> f64 {
        self.obj_int_scale
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::isIntegral (lines 66-66)
    /// C++: `bool isIntegral() const`
    pub fn is_integral(&self) -> bool {
        self.obj_int_scale != 0.0
    }

    /// HiGHS: HighsObjectiveFunction.h::HighsObjectiveFunction::isZero (lines 68-68)
    /// C++: `bool isZero() const`
    pub fn is_zero(&self) -> bool {
        self.objective_nonzeros.is_empty()
    }
}
