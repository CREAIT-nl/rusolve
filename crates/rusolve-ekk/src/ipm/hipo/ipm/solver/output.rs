//! HiGHS source: ipm/hipo/ipm/Solver.cpp (lines 1148-1342)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Interrupt check, log output, solution getters and status predicates.
//! The log lines are HiGHS output and must match byte for byte with
//! `timeless_log` on.

use super::Solver;
use crate::ipm::hipo::Int;
use crate::ipm::hipo::auxiliary::logger::{fix, integer, sci, textline};
use crate::ipm::hipo::ipm::info::Info;
use crate::ipm::hipo::ipm::parameters::{kMaxCorrectors, kMaxIterRefine};
use crate::ipm::hipo::ipm::status::{
    kStatusFailed, kStatusImprecise, kStatusNoProgress, kStatusPDFeas, kStatusSolved, statusString,
};
use crate::lp_data::hconst::{kHighsOffString, kHighsOnString};
use crate::util::highs_utils::{cpp_max, cpp_min};

impl Solver {
    /// HiGHS: Solver.cpp::Solver::checkInterrupt (lines 1148-1156)
    pub fn checkInterrupt(&mut self) -> bool {
        let mut terminate = false;
        let status = self.control_.interruptCheck(self.iter_);
        if status != 0 {
            self.info_.status = status;
            terminate = true;
        }
        terminate
    }

    /// HiGHS: Solver.cpp::Solver::printHeader (lines 1158-1171)
    pub fn printHeader(&self) {
        if self.iter_ % 20 == 0 {
            self.logger_
                .print(" iter       primal obj         dual obj       pinf       dinf       gap");
            if !self.options_.timeless_log {
                self.logger_.print("    time");
            }
            if self.logger_.debug(1) {
                self.logger_.print(
                    "     alpha p/d   sigma af/co   cor  solv    static reg p/d     minT     maxT  (xj * zj / mu)_range_&_num   max_res",
                );
            }
            self.logger_.print("\n");
        }
    }

    /// HiGHS: Solver.cpp::Solver::printOutput (lines 1173-1191)
    pub fn printOutput(&self) {
        self.printHeader();

        let mut line = String::new();
        line.push_str(&integer(self.iter_, 5));
        line.push(' ');
        line.push_str(&sci(self.it_.as_ref().unwrap().pobj, 16, 8));
        line.push(' ');
        line.push_str(&sci(self.it_.as_ref().unwrap().dobj, 16, 8));
        line.push(' ');
        line.push_str(&sci(self.it_.as_ref().unwrap().pinf, 10, 2));
        line.push(' ');
        line.push_str(&sci(self.it_.as_ref().unwrap().dinf, 10, 2));
        line.push(' ');
        line.push_str(&sci(self.it_.as_ref().unwrap().pdgap, 9, 2));
        self.logger_.print(&line);

        if !self.options_.timeless_log {
            let mut line = String::new();
            line.push(' ');
            line.push_str(&fix(self.control_.elapsed(), 7, 1));
            self.logger_.print(&line);
        }

        if self.logger_.debug(1) {
            // PORT NOTE: C++ `it_->data.back()` needs &mut; read the last record via .last().
            let data = self.it_.as_ref().unwrap().data.record.last().unwrap();
            let mut line = String::new();
            line.push(' ');
            line.push_str(&fix(self.alpha_primal_, 6, 2));
            line.push(' ');
            line.push_str(&fix(self.alpha_dual_, 6, 2));
            line.push(' ');
            line.push_str(&fix(data.sigma_aff, 6, 2));
            line.push(' ');
            line.push_str(&fix(data.sigma, 6, 2));
            line.push(' ');
            line.push_str(&integer(data.correctors, 5));
            line.push(' ');
            line.push_str(&integer(data.num_solves, 5));
            line.push(' ');
            line.push_str(&sci(self.regul_.primal, 8, 1));
            line.push(' ');
            line.push_str(&sci(self.regul_.dual, 8, 1));
            line.push(' ');
            line.push_str(&sci(data.min_theta, 8, 1));
            line.push(' ');
            line.push_str(&sci(data.max_theta, 8, 1));
            line.push(' ');
            line.push_str(&sci(data.min_prod, 8, 1));
            line.push(' ');
            line.push_str(&sci(data.max_prod, 8, 1));
            line.push(' ');
            line.push_str(&integer(data.num_small_prod, 4));
            line.push(' ');
            line.push_str(&integer(data.num_large_prod, 4));
            line.push(' ');
            line.push_str(&sci(data.omega, 9, 1));
            self.logger_.print(&line);
        }

        self.logger_.print("\n");
    }

    /// HiGHS: Solver.cpp::Solver::printInfo (lines 1193-1208)
    pub fn printInfo(&self) {
        let mut log_stream = String::new();
        log_stream.push_str("\nRunning HiPO\n");

        // Print number of threads
        if self.options_.parallel == kHighsOffString {
            log_stream.push_str(&textline("Threads:"));
            log_stream.push_str(&integer(1, 0));
            log_stream.push('\n');
        } else {
            // C++ line 1201: highs::parallel::num_threads()
            log_stream.push_str(&textline("Threads:"));
            log_stream.push_str(&integer(
                crate::parallel::highs_parallel::num_threads() as Int,
                0,
            ));
            log_stream.push('\n');
        }

        self.logger_.print(&log_stream);

        // print information about model
        self.model_.print(&self.logger_);
    }

    /// HiGHS: Solver.cpp::Solver::printSummary (lines 1210-1257)
    pub fn printSummary(&self) {
        let mut log_stream = String::new();

        log_stream.push_str("\nSummary\n");
        log_stream.push_str(&crate::ipm::hipo::factorhighs::call_and_time_blas::prof_report());
        if !self.options_.timeless_log {
            log_stream.push_str(&textline("HiPO runtime:"));
            log_stream.push_str(&fix(self.control_.elapsed() - self.start_time_, 0, 2));
            log_stream.push('\n');
        }

        log_stream.push_str(&textline("Status:"));
        log_stream.push_str(statusString(self.info_.status));
        log_stream.push('\n');
        log_stream.push_str(&textline("HiPO iterations:"));
        log_stream.push_str(&integer(self.iter_, 0));
        log_stream.push('\n');
        if self.info_.ipx_used {
            log_stream.push_str(&textline("IPX iterations:"));
            log_stream.push_str(&integer(self.info_.ipx_info.iter, 0));
            log_stream.push('\n');
        }

        if self.logger_.debug(1) {
            log_stream.push_str(&textline("Correctors:"));
            log_stream.push_str(&integer(self.info_.correctors, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Factorisations:"));
            log_stream.push_str(&integer(self.info_.factor_number, 0));
            log_stream.push('\n');
            log_stream.push_str(&textline("Solves:"));
            log_stream.push_str(&integer(self.info_.solve_number, 0));
            log_stream.push('\n');

            if !self.options_.timeless_log {
                log_stream.push_str(&textline("Analyse AS time:"));
                log_stream.push_str(&fix(self.info_.analyse_AS_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("Analyse NE time:"));
                log_stream.push_str(&fix(self.info_.analyse_NE_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("Matrix time:"));
                log_stream.push_str(&fix(self.info_.matrix_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("AS structure time:"));
                log_stream.push_str(&fix(self.info_.AS_structure_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("NE structure time:"));
                log_stream.push_str(&fix(self.info_.NE_structure_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("Factorisation time:"));
                log_stream.push_str(&fix(self.info_.factor_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("Solve time:"));
                log_stream.push_str(&fix(self.info_.solve_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("Residual time:"));
                log_stream.push_str(&fix(self.info_.residual_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("Omega time:"));
                log_stream.push_str(&fix(self.info_.omega_time, 0, 2));
                log_stream.push('\n');
                log_stream.push_str(&textline("Avg time per factorisation:"));
                log_stream.push_str(&sci(
                    self.info_.factor_time / self.info_.factor_number as f64,
                    0,
                    2,
                ));
                log_stream.push('\n');
                log_stream.push_str(&textline("Avg time per solve:"));
                log_stream.push_str(&sci(
                    self.info_.solve_time / self.info_.solve_number as f64,
                    0,
                    2,
                ));
                log_stream.push('\n');
            }
        }

        self.logger_.print(&log_stream);
    }

    /// HiGHS: Solver.cpp::Solver::getOriginalDims (lines 1259-1262)
    pub fn getOriginalDims(&self, num_row: &mut Int, num_col: &mut Int) {
        *num_row = self.model_.m_orig();
        *num_col = self.model_.n_orig();
    }

    /// HiGHS: Solver.cpp::Solver::getInfo (lines 1263-1263)
    pub fn getInfo(&self) -> &Info {
        &self.info_
    }

    /// HiGHS: Solver.cpp::Solver::getInteriorSolution (lines 1264-1276)
    pub fn getInteriorSolution(
        &self,
        x: &mut Vec<f64>,
        xl: &mut Vec<f64>,
        xu: &mut Vec<f64>,
        slack: &mut Vec<f64>,
        y: &mut Vec<f64>,
        zl: &mut Vec<f64>,
        zu: &mut Vec<f64>,
    ) {
        // prepare and return solution with internal format

        if !self.info_.ipx_used {
            self.model_
                .postprocess(x, xl, xu, slack, y, zl, zu, self.it_.as_ref().unwrap());
        } else {
            self.ipx_lps_.get_interior_solution(
                Some(x.as_mut_slice()),
                Some(xl.as_mut_slice()),
                Some(xu.as_mut_slice()),
                Some(slack.as_mut_slice()),
                Some(y.as_mut_slice()),
                Some(zl.as_mut_slice()),
                Some(zu.as_mut_slice()),
            );
        }
    }

    /// HiGHS: Solver.cpp::Solver::getBasicSolution (lines 1278-1284)
    pub fn getBasicSolution(
        &self,
        x: &mut Vec<f64>,
        slack: &mut Vec<f64>,
        y: &mut Vec<f64>,
        z: &mut Vec<f64>,
        cbasis: &mut [Int],
        vbasis: &mut [Int],
    ) -> Int {
        // interface to ipx getBasicSolution
        self.ipx_lps_.get_basic_solution(
            Some(x.as_mut_slice()),
            Some(slack.as_mut_slice()),
            Some(y.as_mut_slice()),
            Some(z.as_mut_slice()),
            Some(cbasis),
            Some(vbasis),
        )
    }

    /// HiGHS: Solver.cpp::Solver::maxCorrectors (lines 1286-1321)
    pub fn maxCorrectors(&mut self) {
        if kMaxCorrectors > 0 {
            // Compute estimate of effort to factorise and solve

            // Effort to factorise depends on the number of flops
            let fact_effort =
                self.LS_.as_ref().unwrap().flops() + 100.0 * self.LS_.as_ref().unwrap().spops();

            // Effort to solve depends on the number of nonzeros of L multiplied by 2,
            // because there are two sweeps through L (forward and backward).
            let solv_effort = 2.0 * self.LS_.as_ref().unwrap().nz();

            // The factorise phase uses BLAS-3 and can be parallelised, the solve phase
            // uses BLAS-2 and cannot be parallelised. To account for this, the
            // factorisation effort is multiplied by a coefficient < 1, estimated
            // empirically.
            let alpha = 1.0 / 112.0;

            let ratio = alpha * fact_effort / solv_effort;

            // At each ipm iteration, there are up to (1+k) directions computed, where k
            // is the number of correctors. Each direction requires up (1+f) solves,
            // where f is the number of refinement steps. So, up to (1+k)(1+f) solves
            // are performed per iteration. However, not all refinement steps are used
            // all the time, so use f/2.
            // Therefore, we want (1+k)(1+f/2) < ratio.

            let thresh = ratio / (1.0 + kMaxIterRefine as f64 / 2.0) - 1.0;

            self.info_.correctors = thresh.floor() as Int;
            self.info_.correctors = cpp_max(self.info_.correctors, 1);
            self.info_.correctors = cpp_min(self.info_.correctors, kMaxCorrectors);
        } else {
            self.info_.correctors = -kMaxCorrectors;
        }
    }

    /// HiGHS: Solver.cpp::Solver::statusIsSolved (lines 1323-1323)
    pub fn statusIsSolved(&self) -> bool {
        self.info_.status >= kStatusSolved
    }

    /// HiGHS: Solver.cpp::Solver::statusIsStopped (lines 1324-1324)
    pub fn statusIsStopped(&self) -> bool {
        self.info_.status < kStatusFailed
    }

    /// HiGHS: Solver.cpp::Solver::statusIsFailed (lines 1325-1327)
    pub fn statusIsFailed(&self) -> bool {
        self.info_.status >= kStatusFailed && self.info_.status < kStatusSolved
    }

    /// HiGHS: Solver.cpp::Solver::statusAllowsCrossover (lines 1328-1330)
    pub fn statusAllowsCrossover(&self) -> bool {
        self.info_.status >= kStatusPDFeas
    }

    /// HiGHS: Solver.cpp::Solver::statusNeedsRefinement (lines 1331-1333)
    pub fn statusNeedsRefinement(&self) -> bool {
        self.info_.status == kStatusNoProgress || self.info_.status == kStatusImprecise
    }

    /// HiGHS: Solver.cpp::Solver::refinementIsOn (lines 1334-1336)
    pub fn refinementIsOn(&self) -> bool {
        self.options_.refine_with_ipx && !self.model_.qp()
    }

    /// HiGHS: Solver.cpp::Solver::crossoverIsOn (lines 1337-1339)
    pub fn crossoverIsOn(&self) -> bool {
        self.options_.crossover == kHighsOnString && !self.model_.qp()
    }

    /// HiGHS: Solver.cpp::Solver::solved (lines 1340-1340)
    pub fn solved(&self) -> bool {
        self.statusIsSolved()
    }

    /// HiGHS: Solver.cpp::Solver::stopped (lines 1341-1341)
    pub fn stopped(&self) -> bool {
        self.statusIsStopped()
    }

    /// HiGHS: Solver.cpp::Solver::failed (lines 1342-1342)
    pub fn failed(&self) -> bool {
        self.statusIsFailed()
    }
}
