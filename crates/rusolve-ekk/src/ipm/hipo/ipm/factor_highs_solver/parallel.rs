//! HiGHS source: ipm/hipo/ipm/FactorHiGHSSolver.cpp (lines 202-355, 526-659)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! `chooseNla` (`hipo_system=choose`: analyse NE, then AS, compare),
//! `setNla` (explicit `hipo_system`) and `setParallel` (tree / node
//! parallelism; with one thread it resolves to serial and only logs).
//! With `options.parallel == "off"` the two steps of `chooseNla` run in
//! sequence (AS analyse, then NE structure); otherwise they are spawned in a
//! `TaskGroup` (lines 275-280) on disjoint borrows of `kkt_` and `info_`.
//! `chooseNla` sits here with its caller `setNla` (not in `choose.rs`).
//! `usingAppleBlas` is complete: the BLAS is [`crate::ipm::hipo::blas`]
//! (refblas), never Apple Accelerate.

use super::FactorHiGHSSolver;
use crate::ipm::hipo::Int;
use crate::ipm::hipo::Int64;
use crate::ipm::hipo::auxiliary::logger::{Logger, textline};
use crate::ipm::hipo::factorhighs::symbolic::Symbolic;
use crate::ipm::hipo::ipm::kkt_matrix::{KktAsPart, KktNePart};
use crate::ipm::hipo::ipm::linear_solver::LsCtx;
use crate::ipm::hipo::ipm::model::Model;
use crate::ipm::hipo::ipm::options::Options;
use crate::ipm::hipo::ipm::parameters::*;
use crate::ipm::hipo::ipm::status::*;
use crate::lp_data::hconst::{kHighsChooseString, kHighsIInf, kHighsOffString, kHighsOnString};
use crate::lp_data::highs_options::{
    kHipoAugmentedString, kHipoBothString, kHipoNodeString, kHipoNormalEqString, kHipoTreeString,
};
use std::sync::atomic::Ordering;

/// Ours (no C++ counterpart): the `const Model& model_` that both `chooseNla` tasks read.
/// `Model` is not `Sync` only because its `Preprocessor` holds `Box<dyn PreprocessAction>`,
/// which neither task touches.
#[derive(Clone, Copy)]
struct ModelRef<'a>(&'a Model);
// SAFETY: C++ `const Model&` shared by the AS and NE tasks of `chooseNla` (lines 275-280):
// both only read it (A, Q, qp, m/n); nothing writes the model while the group runs.
unsafe impl Sync for ModelRef<'_> {}
// SAFETY: as above; a copy of the shared reference is moved into a task.
unsafe impl Send for ModelRef<'_> {}
impl<'a> ModelRef<'a> {
    // A method call, so that closures capture the whole `ModelRef`, not its field.
    fn get(self) -> &'a Model {
        self.0
    }
}

impl FactorHiGHSSolver {
    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::chooseNla (lines 202-355)
    pub fn chooseNla(&mut self, ctx: &mut LsCtx<'_>) -> Int {
        // Choose whether to use augmented system or normal equations.

        let mut symb_NE = Symbolic::new();
        let mut symb_AS = Symbolic::new();
        let mut failure_NE = false;
        let mut failure_AS = false;
        let mut overflow_NE = false;
        let mut overflow_AS = false;

        let AS_possible = true;
        let NE_possible = !(ctx.model.nonSeparableQp() || ctx.model.m() == 0);

        // C++: `Int64 > Int64 * double`, both sides converted to double.
        let expect_AS_much_cheaper =
            (ctx.model.nzNElb() as f64) > (ctx.model.nzAS() as f64) * kNzBoundsRatio;
        let expect_NE_much_cheaper =
            (ctx.model.nzAS() as f64) > (ctx.model.nzNEub() as f64) * kNzBoundsRatio;

        let skip_AS = NE_possible && expect_NE_much_cheaper;
        let skip_NE = AS_possible && expect_AS_much_cheaper;

        let mut allow_skip_AS = true;
        let mut allow_skip_NE = true;

        // The C++ lambdas capture `this` (`kkt_`, `info_`, `logger_`) by reference; the two
        // tasks touch disjoint members of `kkt_` and `info_` (AS: ptrAS/rowsAS/valAS,
        // AS_structure_time, analyse_AS_time, FH_, ordering_AS_; NE: ptrNE/rowsNE/ptrA_rw/
        // idxA_rw/corr_A, NE_structure_time) and share only `model_`, `options_`,
        // `logger_` (read) and the atomic `NE_nz_limit`. Here those members are explicit
        // arguments, borrowed apart.
        let model = ModelRef(ctx.model);
        let logger = ctx.logger;

        let run_structure_NE = |kkt_NE: &mut KktNePart<'_>,
                                NE_structure_time: &mut f64,
                                failure_NE: &mut bool,
                                overflow_NE: &mut bool,
                                allow_skip_NE: bool| {
            if !NE_possible || (skip_NE && allow_skip_NE) {
                *failure_NE = true;
                logger.printInfo("NE skipped\n");
            } else {
                let status = kkt_NE.buildNEstructure(model.get(), NE_structure_time, logger);
                if status != 0 {
                    *failure_NE = true;
                    if status == kStatusOverflow {
                        logger.printInfo("Integer overflow forming NE matrix\n");
                        *overflow_NE = true;
                    }
                }
            }
        };

        let run_analyse_AS = |this: &mut FactorHiGHSSolver,
                              options: &Options,
                              kkt_AS: &mut KktAsPart<'_>,
                              AS_structure_time: &mut f64,
                              analyse_AS_time: &mut f64,
                              symb_AS: &mut Symbolic,
                              failure_AS: &mut bool,
                              overflow_AS: &mut bool,
                              allow_skip_AS: bool| {
            if !AS_possible || (skip_AS && allow_skip_AS) {
                *failure_AS = true;
                logger.printInfo("AS skipped\n");
            } else {
                let mut AS_status = kkt_AS.buildASstructure(model.get(), AS_structure_time, logger);
                if AS_status == 0 {
                    AS_status = this.analyseAS_part(
                        kkt_AS,
                        model.get(),
                        options,
                        logger,
                        analyse_AS_time,
                        symb_AS,
                    );
                }
                if AS_status != 0 {
                    *failure_AS = true;
                }
                if AS_status == kStatusOverflow {
                    logger.printInfo("Integer overflow forming AS matrix\n");
                    *overflow_AS = true;
                }

                // If NE has more nonzeros than the factor of AS, then it's likely that AS
                // will be preferred, so stop computation of NE.
                let mut NE_nz_limit = (symb_AS.nz() as f64 * kSymbNzMult) as Int64;
                if *failure_AS || NE_nz_limit > kHighsIInf as Int64 {
                    NE_nz_limit = kHighsIInf as Int64;
                }
                kkt_AS.NE_nz_limit.store(NE_nz_limit, Ordering::Relaxed);
            }
        };

        // rusolve (ipm_hipo_amd_first): analyse NE first. When its
        // predicted factor is cheap (flops <= 5e10) AS is not analysed at all: on rail4284 the
        // AS analyse (dimension 1.1M, AMD plus Metis) took 95 s against 0.3 s for NE.
        let mut ne_done = false;
        let mut as_skipped_by_ne = false;
        // NE upper bound capped: without AS's NE_nz_limit a huge NE is built and ordered
        // (Linf_520c was killed at 8 GB)
        if ctx.options.amd_first
            && NE_possible
            && !skip_NE
            && !skip_AS
            && (ctx.model.nzNEub() as f64) <= 1.5e7
        {
            {
                let (_, mut kkt_NE) = ctx.kkt.split_AS_NE();
                run_structure_NE(
                    &mut kkt_NE,
                    &mut ctx.info.NE_structure_time,
                    &mut failure_NE,
                    &mut overflow_NE,
                    allow_skip_NE,
                );
            }
            if !failure_NE && self.analyseNE(ctx, &mut symb_NE) != 0 {
                failure_NE = true;
            }
            ne_done = true;
            if !failure_NE && symb_NE.flops() <= 5.0e10 {
                failure_AS = true;
                as_skipped_by_ne = true;
                logger.printInfo("AS skipped (NE analysed first, cheap)\n");
            }
        }

        // In parallel, run AS analyse and build NE structure. NE analyse runs only
        // after AS analyse is finished, so that it can be skipped based on the number
        // of nz of NE matrix and AS factor.
        if !as_skipped_by_ne {
            let options: &Options = ctx.options;
            let (mut kkt_AS, mut kkt_NE) = ctx.kkt.split_AS_NE();
            let info = &mut *ctx.info;
            let AS_structure_time = &mut info.AS_structure_time;
            let analyse_AS_time = &mut info.analyse_AS_time;
            let NE_structure_time = &mut info.NE_structure_time;
            if options.parallel == kHighsOffString {
                run_analyse_AS(
                    self,
                    options,
                    &mut kkt_AS,
                    AS_structure_time,
                    analyse_AS_time,
                    &mut symb_AS,
                    &mut failure_AS,
                    &mut overflow_AS,
                    allow_skip_AS,
                );
                if !ne_done {
                    run_structure_NE(
                        &mut kkt_NE,
                        NE_structure_time,
                        &mut failure_NE,
                        &mut overflow_NE,
                        allow_skip_NE,
                    );
                }
            } else {
                // `taskWait` runs the spawned tasks LIFO on this worker: at one worker
                // NE structure runs before AS analyse, with `NE_nz_limit` still
                // kHighsIInf.
                crate::parallel::highs_parallel::task_group(|tg| {
                    tg.spawn(|| {
                        run_analyse_AS(
                            self,
                            options,
                            &mut kkt_AS,
                            AS_structure_time,
                            analyse_AS_time,
                            &mut symb_AS,
                            &mut failure_AS,
                            &mut overflow_AS,
                            allow_skip_AS,
                        )
                    });
                    tg.spawn(|| {
                        if !ne_done {
                            run_structure_NE(
                                &mut kkt_NE,
                                NE_structure_time,
                                &mut failure_NE,
                                &mut overflow_NE,
                                allow_skip_NE,
                            )
                        }
                    });
                    tg.task_wait();
                });
            }
        }

        // if NE was skipped but AS failed, use NE
        if skip_NE && failure_AS {
            allow_skip_NE = false;
            let (_, mut kkt_NE) = ctx.kkt.split_AS_NE();
            run_structure_NE(
                &mut kkt_NE,
                &mut ctx.info.NE_structure_time,
                &mut failure_NE,
                &mut overflow_NE,
                allow_skip_NE,
            );
        }

        // run_analyse_NE (lines 243-247)
        if !failure_NE && !ne_done {
            let NE_status = self.analyseNE(ctx, &mut symb_NE);
            if NE_status != 0 {
                failure_NE = true;
            }
        }

        // if AS was skipped but NE failed, use AS
        if skip_AS && failure_NE {
            allow_skip_AS = false;
            let (mut kkt_AS, _) = ctx.kkt.split_AS_NE();
            let info = &mut *ctx.info;
            run_analyse_AS(
                self,
                ctx.options,
                &mut kkt_AS,
                &mut info.AS_structure_time,
                &mut info.analyse_AS_time,
                &mut symb_AS,
                &mut failure_AS,
                &mut overflow_AS,
                allow_skip_AS,
            );
        }

        let mut status = kStatusOk;

        let mut log_stream = String::new();

        // Decision may be forced by failures
        if failure_NE && !failure_AS {
            ctx.options.nla = kHipoAugmentedString.to_string();
            log_stream.push_str(&textline("Newton system:"));
            log_stream.push_str("AS preferred (NE failed)\n");
        } else if failure_AS && !failure_NE {
            ctx.options.nla = kHipoNormalEqString.to_string();
            log_stream.push_str(&textline("Newton system:"));
            log_stream.push_str("NE preferred (AS failed)\n");
        } else if failure_AS && failure_NE {
            if overflow_AS && overflow_NE {
                status = kStatusOverflow;
            } else {
                status = kStatusErrorAnalyse;
            }

            ctx.logger.printe("Both NE and AS failed analyse phase\n");
        } else {
            // Total number of operations, given by dense flops and sparse indexing
            // operations, weighted with an empirical factor
            let ops_NE = symb_NE.flops() + symb_NE.spops() * kSpopsWeight;
            let ops_AS = symb_AS.flops() + symb_AS.spops() * kSpopsWeight;

            // Average size of supernodes
            let sn_size_NE = (symb_NE.size() as f64) / (symb_NE.sn() as f64);
            let sn_size_AS = (symb_AS.size() as f64) / (symb_AS.sn() as f64);

            let ratio_ops = ops_NE / ops_AS;
            let ratio_sn = sn_size_AS / sn_size_NE;

            let NE_much_more_expensive = ratio_ops > kRatioOpsThresh;
            let AS_not_too_expensive = ratio_ops > 1.0 / kRatioOpsThresh;
            let sn_AS_larger_than_NE = ratio_sn > kRatioSnThresh;

            if NE_much_more_expensive || (sn_AS_larger_than_NE && AS_not_too_expensive) {
                ctx.options.nla = kHipoAugmentedString.to_string();
                log_stream.push_str(&textline("Newton system:"));
                log_stream.push_str("AS preferred\n");
            } else {
                ctx.options.nla = kHipoNormalEqString.to_string();
                log_stream.push_str(&textline("Newton system:"));
                log_stream.push_str("NE preferred\n");
            }
        }

        ctx.logger.print(&log_stream);

        if status == kStatusOk {
            if ctx.options.nla == kHipoAugmentedString {
                self.S_ = symb_AS;
                ctx.kkt.freeNEmemory();
            } else {
                self.S_ = symb_NE;
                ctx.kkt.freeASmemory();
            }
        }

        status
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::setNla (lines 526-574)
    pub fn setNla(&mut self, ctx: &mut LsCtx<'_>) -> Int {
        let mut log_stream = String::new();

        if ctx.options.nla == kHipoNormalEqString && ctx.model.nonSeparableQp() {
            ctx.logger
                .printw("Normal equations not available for non-separable QP\n");
            ctx.options.nla = kHighsChooseString.to_string();
        }

        if ctx.options.nla == kHipoAugmentedString {
            let mut status = ctx.kkt.buildASstructure(ctx.model, ctx.info, ctx.logger);
            if status == 0 {
                // C++ passes S_ by reference; analyseAS never reads self.S_, so the
                // value is taken out for the call and put back (borrow split).
                let mut S_tmp = std::mem::take(&mut self.S_);
                status = self.analyseAS(ctx, &mut S_tmp);
                self.S_ = S_tmp;
            }
            if status == kStatusOverflow {
                ctx.logger.printe("AS requested, integer overflow\n");
                return kStatusOverflow;
            } else if status != 0 {
                ctx.logger.printe("AS requested, failed analyse phase\n");
                return kStatusErrorAnalyse;
            }
            log_stream.push_str(&textline("Newton system:"));
            log_stream.push_str("AS requested\n");
        } else if ctx.options.nla == kHipoNormalEqString {
            let mut status = ctx.kkt.buildNEstructure(ctx.model, ctx.info, ctx.logger);
            if status == 0 {
                // Same borrow split as analyseAS above; analyseNE never reads self.S_.
                let mut S_tmp = std::mem::take(&mut self.S_);
                status = self.analyseNE(ctx, &mut S_tmp);
                self.S_ = S_tmp;
            }
            if status == kStatusOverflow {
                ctx.logger.printe("NE requested, integer overflow\n");
                return kStatusOverflow;
            } else if status != 0 {
                ctx.logger.printe("NE requested, failed analyse phase\n");
                return kStatusErrorAnalyse;
            }
            log_stream.push_str(&textline("Newton system:"));
            log_stream.push_str("NE requested\n");
        } else if ctx.options.nla == kHighsChooseString {
            let status = self.chooseNla(ctx);
            if status != 0 {
                return status;
            }
        } else {
            debug_assert!(false);
        }

        if ctx.logger.debug(1) {
            log_stream.push_str(&textline("Ordering:"));
            if ctx.options.nla == kHipoAugmentedString {
                log_stream.push_str(&self.ordering_AS_);
            } else {
                log_stream.push_str(&self.ordering_NE_);
            }
            log_stream.push('\n');
        }
        ctx.logger.print(&log_stream);

        ctx.kkt.iperm = self.S_.iperm().to_vec();

        kStatusOk
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::setParallel (lines 580-659)
    pub fn setParallel(&mut self, ctx: &mut LsCtx<'_>) {
        // Set parallel options
        let mut parallel_tree = false;
        let mut parallel_node = false;

        let mut log_stream = String::new();
        log_stream.push_str(&textline("Parallelism:"));

        if ctx.options.parallel == kHighsOffString {
            log_stream.push_str("None requested\n");
        } else if ctx.options.parallel == kHighsOnString {
            if ctx.options.parallel_type == kHipoBothString {
                parallel_tree = true;
                parallel_node = true;
                log_stream.push_str("Full requested\n");
            } else if ctx.options.parallel_type == kHipoTreeString {
                parallel_tree = true;
                log_stream.push_str("Tree requested\n");
            } else if ctx.options.parallel_type == kHipoNodeString {
                parallel_node = true;
                log_stream.push_str("Node requested\n");
            } else {
                debug_assert!(false);
            }
        } else if ctx.options.parallel == kHighsChooseString {
            if usingAppleBlas() {
                // Blas on Apple do not work well with parallel_node, but parallel_tree
                // seems to always be beneficial.
                parallel_node = false;
                parallel_tree = true;
            } else {
                // Otherwise, parallel_node is active because it is triggered only if the
                // frontal matrix is large enough anyway.
                parallel_node = true;

                // parallel_tree instead is chosen with a heuristic

                let tree_speedup = self.S_.flops() / self.S_.critops();
                let sn_size = (self.S_.size() as f64) / (self.S_.sn() as f64);

                let enough_sn = self.S_.sn() > kMinNumberSn;
                let enough_flops = self.S_.flops() > kLargeFlopsThresh;
                let speedup_is_large = tree_speedup > kLargeSpeedupThresh;
                let sn_are_large = sn_size > kLargeSnThresh;
                let sn_are_not_small = sn_size > kSmallSnThresh;

                // parallel_tree is active if the supernodes are large, or if there is a
                // large expected speedup and the supernodes are not too small, provided
                // that the number of flops and supernodes is not too small.
                if enough_sn
                    && enough_flops
                    && (sn_are_large || (speedup_is_large && sn_are_not_small))
                {
                    parallel_tree = true;
                }
            }

            // If serial memory is too large, switch off tree parallelism to avoid
            // running out of memory
            let num_GB = self.S_.storage() / 1024.0 / 1024.0 / 1024.0;
            if num_GB > kLargeStorageGB {
                parallel_tree = false;
            }

            // switch off tree parallelism if depth of recursion is too large
            if (self.S_.depth() as f64) > kMaxTreeDepth {
                parallel_tree = false;
            }

            if parallel_tree && parallel_node {
                log_stream.push_str("Full preferred\n");
            } else if parallel_tree && !parallel_node {
                log_stream.push_str("Tree preferred\n");
            } else if !parallel_tree && parallel_node {
                log_stream.push_str("Node preferred\n");
            } else {
                log_stream.push_str("None preferred\n");
            }
        } else {
            debug_assert!(false);
        }

        ctx.logger.print(&log_stream);
        self.S_.setParallel(parallel_tree, parallel_node);
    }
}

/// HiGHS: FactorHiGHSSolver.cpp::usingAppleBlas (lines 576-578)
pub fn usingAppleBlas() -> bool {
    false
}

#[cfg(test)]
mod tests {
    //! Unit tests, hand-derived from FactorHiGHSSolver.cpp
    //! `setParallel` (lines 580-659).
    use super::*;
    use crate::ipm::hipo::ipm::info::Info;
    use crate::ipm::hipo::ipm::ipm_data::IpmData;
    use crate::ipm::hipo::ipm::kkt_matrix::KktMatrix;
    use crate::ipm::hipo::ipm::model::Model;
    use crate::ipm::hipo::ipm::options::Options;
    use crate::ipm::hipo::ipm::parameters::Regularisation;

    /// Run `setParallel` on a solver whose `S_` is `S` and return the
    /// resulting (parallel_tree, parallel_node).
    fn run(parallel: &str, parallel_type: &str, S: Symbolic) -> (bool, bool) {
        let mut kkt = KktMatrix::new();
        let mut options = Options::default();
        options.parallel = parallel.to_string();
        options.parallel_type = parallel_type.to_string();
        let model = Model::default();
        let regul = Regularisation::default();
        let mut info = Info::default();
        let mut data = IpmData::default();
        let logger = Logger::default();
        let mut s = FactorHiGHSSolver::new(&options, &logger);
        s.S_ = S;
        let mut ctx = LsCtx {
            kkt: &mut kkt,
            options: &mut options,
            model: &model,
            regul: &regul,
            info: &mut info,
            data: &mut data,
            logger: &logger,
        };
        s.setParallel(&mut ctx);
        (s.S_.parTree(), s.S_.parNode())
    }

    /// Large supernodes, many flops: the `choose` heuristic picks tree.
    fn big_symbolic() -> Symbolic {
        let mut S = Symbolic::new();
        S.n_ = 1000;
        S.sn_ = 20; // sn_size 50 > kLargeSnThresh, sn 20 > kMinNumberSn
        S.flops_ = 1e8; // > kLargeFlopsThresh
        S.critops_ = 1e7;
        S.serial_storage_ = 0.0;
        S.tree_depth_ = 5;
        S
    }

    #[test]
    fn requested_modes() {
        let S = big_symbolic;
        assert_eq!(run(kHighsOffString, kHipoBothString, S()), (false, false));
        assert_eq!(run(kHighsOnString, kHipoBothString, S()), (true, true));
        assert_eq!(run(kHighsOnString, kHipoTreeString, S()), (true, false));
        assert_eq!(run(kHighsOnString, kHipoNodeString, S()), (false, true));
    }

    #[test]
    fn choose_heuristic() {
        // Empty symbolic: sn 0 is not > kMinNumberSn, so node only.
        assert_eq!(
            run(kHighsChooseString, kHipoBothString, Symbolic::new()),
            (false, true)
        );
        assert_eq!(
            run(kHighsChooseString, kHipoBothString, big_symbolic()),
            (true, true)
        );
        // Medium supernodes (10): need speedup > kLargeSpeedupThresh.
        let mut S = big_symbolic();
        S.n_ = 200;
        S.critops_ = 1e8; // speedup 1.0, not > 1.0
        assert_eq!(
            run(kHighsChooseString, kHipoBothString, S.clone()),
            (false, true)
        );
        S.critops_ = 5e7; // speedup 2.0
        assert_eq!(run(kHighsChooseString, kHipoBothString, S), (true, true));
        // Storage above kLargeStorageGB switches tree off (line 637-640).
        let mut S = big_symbolic();
        S.serial_storage_ = 21.0 * 1024.0 * 1024.0 * 1024.0;
        assert_eq!(run(kHighsChooseString, kHipoBothString, S), (false, true));
        // Depth above kMaxTreeDepth switches tree off (line 643).
        let mut S = big_symbolic();
        S.tree_depth_ = 1001;
        assert_eq!(
            run(kHighsChooseString, kHipoBothString, S.clone()),
            (false, true)
        );
        S.tree_depth_ = 1000;
        assert_eq!(run(kHighsChooseString, kHipoBothString, S), (true, true));
    }

    // Unit tests, hand-derived from FactorHiGHSSolver.cpp `chooseNla`
    // (lines 202-355), KktMatrix.cpp 170-173 and FactorHiGHSSolver.cpp 381.

    use crate::io::highs_io::{HighsLogOptions, kHighsLogDevLevelInfo};
    use crate::lp_data::hconst::kHighsInf;
    use crate::lp_data::highs_callback::{
        HighsCallbackFunctionType, HighsCallbackInput, HighsCallbackOutput,
    };
    use crate::lp_data::highs_lp::HighsLp;
    use crate::model::highs_hessian::HighsHessian;
    use crate::parallel::highs_parallel::initialize_scheduler;
    use crate::parallel::highs_task_executor::HighsTaskExecutor;
    use std::sync::{Arc, Mutex};

    /// Equality LP with `m` rows: an identity block plus, when `dense`, one column with
    /// an entry in every row (so A*A^T is dense while AS stays an arrow).
    fn lp(m: Int, dense: bool) -> HighsLp {
        let n = if dense { m + 1 } else { m };
        let mut lp = HighsLp {
            num_col_: n,
            num_row_: m,
            col_cost_: vec![1.0; n as usize],
            col_lower_: vec![0.0; n as usize],
            col_upper_: vec![kHighsInf; n as usize],
            row_lower_: vec![1.0; m as usize],
            row_upper_: vec![1.0; m as usize],
            ..Default::default()
        };
        let mut start = vec![0];
        let mut index = Vec::new();
        let mut value = Vec::new();
        for i in 0..m {
            index.push(i);
            value.push(1.0);
            if !dense && i + 1 < m {
                // tridiagonal-ish sparse coupling
                index.push(i + 1);
                value.push(0.5);
            }
            start.push(index.len() as Int);
        }
        if dense {
            for i in 0..m {
                index.push(i);
                value.push(2.0);
            }
            start.push(index.len() as Int);
        }
        lp.a_matrix_.num_col_ = n;
        lp.a_matrix_.num_row_ = m;
        lp.a_matrix_.start_ = start;
        lp.a_matrix_.index_ = index;
        lp.a_matrix_.value_ = value;
        lp
    }

    struct Outcome {
        status: Int,
        nla: String,
        NE_nz_limit: Int64,
        S_nz: Int64,
        iperm: Vec<Int>,
        log: String,
    }

    /// Run `chooseNla` with `parallel` on `lp` (or the default model when `None`) on the
    /// current thread.
    fn choose(lp: Option<&HighsLp>, parallel: &str) -> Outcome {
        let mut model = Model::default();
        if let Some(lp) = lp {
            assert_eq!(model.init(lp, &HighsHessian::default()), 0);
        }
        let regul = Regularisation::default();
        let out = Arc::new(Mutex::new(String::new()));
        let sink = out.clone();
        let cb = HighsCallbackFunctionType::new(
            move |_t: Int, msg: &str, _o: &HighsCallbackOutput, _i: &mut HighsCallbackInput| {
                sink.lock().unwrap().push_str(msg);
            },
        );
        let logger = Logger {
            log_options_: Some(HighsLogOptions {
                log_to_console: false,
                log_dev_level: kHighsLogDevLevelInfo,
                user_callback: Some(cb),
                user_callback_active: true,
                ..Default::default()
            }),
        };
        let mut kkt = KktMatrix::new();
        let mut options = Options::default();
        options.nla = kHighsChooseString.to_string();
        options.parallel = parallel.to_string();
        let mut info = Info::default();
        let mut data = IpmData::default();
        let mut s = FactorHiGHSSolver::new(&options, &logger);
        let mut ctx = LsCtx {
            kkt: &mut kkt,
            options: &mut options,
            model: &model,
            regul: &regul,
            info: &mut info,
            data: &mut data,
            logger: &logger,
        };
        let status = s.chooseNla(&mut ctx);
        let nla = ctx.options.nla.clone();
        let NE_nz_limit = ctx.kkt.NE_nz_limit.load(Ordering::Relaxed);
        let log = out.lock().unwrap().clone();
        Outcome {
            status,
            nla,
            NE_nz_limit,
            S_nz: s.S_.nz(),
            iperm: s.S_.iperm().to_vec(),
            log,
        }
    }

    /// `chooseNla` on a fresh thread with its own `threads`-worker scheduler.
    fn choose_on_thread(lp: Option<HighsLp>, parallel: &'static str, threads: i32) -> Outcome {
        std::thread::spawn(move || {
            initialize_scheduler(threads);
            let r = choose(lp.as_ref(), parallel);
            HighsTaskExecutor::shutdown(true);
            r
        })
        .join()
        .unwrap()
    }

    /// The serial branch (lines 271-273).
    fn choose_serial(lp: Option<HighsLp>) -> Outcome {
        choose_on_thread(lp, kHighsOffString, 1)
    }

    /// The TaskGroup branch (lines 274-280) at `threads` workers.
    fn choose_in_group(lp: Option<HighsLp>, threads: i32) -> Outcome {
        choose_on_thread(lp, kHighsChooseString, threads)
    }

    fn assert_same(serial: &Outcome, group: &Outcome) {
        assert_eq!(serial.status, group.status);
        assert_eq!(serial.nla, group.nla);
        assert_eq!(serial.NE_nz_limit, group.NE_nz_limit);
        assert_eq!(serial.S_nz, group.S_nz);
        assert_eq!(serial.iperm, group.iperm);
    }

    /// HiGHS: FactorHiGHSSolver.cpp::FactorHiGHSSolver::chooseNla (lines 202-355): the
    /// TaskGroup branch gives the serial result at every thread count.
    #[test]
    fn w5_hipo_nla_choose_nla_task_group() {
        // m = 0: NE not possible, AS structure has no rows; both fail (line 316).
        // Sparse LP: both analyses succeed and the ratios decide (lines 322-346).
        for lp in [None, Some(lp(30, false))] {
            let serial = choose_serial(lp.clone());
            for threads in [1, 2, 4, 8] {
                for _ in 0..2 {
                    assert_same(&serial, &choose_in_group(lp.clone(), threads));
                }
            }
        }
        let serial = choose_serial(None);
        assert_eq!(serial.status, kStatusErrorAnalyse);
        assert_eq!(serial.NE_nz_limit, kHighsIInf as Int64);
        let serial = choose_serial(Some(lp(30, false)));
        assert_eq!(serial.status, kStatusOk);
    }

    /// NE has more than `symb_AS.nz() * kSymbNzMult` entries: serially AS runs first and
    /// NE fails in `buildNEstructure` (KktMatrix.cpp 170-173, log line 234); in the group at
    /// one worker `taskWait` runs NE structure first (LIFO) with the limit still inf, and NE
    /// fails at `chooseOrdering` line 381 ("NE interrupted"). Same nla, status and `S_`.
    #[test]
    fn w5_hipo_nla_ne_limit_log_line() {
        let dense = lp(200, true);
        let serial = choose_serial(Some(dense.clone()));
        assert_eq!(serial.status, kStatusOk);
        assert_eq!(serial.nla, kHipoAugmentedString);
        assert!(serial.NE_nz_limit < kHighsIInf as Int64);
        assert!(serial.log.contains("Integer overflow forming NE matrix"));
        assert!(!serial.log.contains("NE interrupted"));
        assert!(serial.log.contains("AS preferred (NE failed)"));

        let one = choose_in_group(Some(dense.clone()), 1);
        assert_same(&serial, &one);
        assert!(one.log.contains("NE interrupted"));
        assert!(!one.log.contains("Integer overflow forming NE matrix"));
        assert!(one.log.contains("AS preferred (NE failed)"));
        // Order at one worker: NE structure is built before AS analyse.
        let ne = one.log.find("Building NE structure").unwrap();
        let as_ = one.log.find("Building AS structure").unwrap();
        assert!(ne < as_);

        // At more workers only which of the two log lines appears is racy.
        for threads in [2, 4, 8] {
            for _ in 0..2 {
                assert_same(&serial, &choose_in_group(Some(dense.clone()), threads));
            }
        }
    }
}
