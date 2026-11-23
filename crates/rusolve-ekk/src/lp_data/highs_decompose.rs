//! rusolve: solve the connected components of a MIP as separate MIPs.
//!
//! Not in HiGHS. Rows and columns that share no coefficient form independent subproblems; the
//! optimum is the sum of the component optima and a branch-and-bound tree over the sum is the
//! sum of the component trees rather than their product. The LP of a block-diagonal model is
//! also solved faster per block than as a whole (neos-5052403-cygnet, six blocks: the whole
//! root LP does not finish in 300 s, each block does in 40 to 100 s).
//! Textbook source: Achterberg, "Constraint Integer Programming" (2007) section 10.1 on
//! decomposition of independent subproblems; Gamrath et al., "Exploiting solving phases"
//! (2020) for SCIP's component presolver (`presol_components`, Apache-2.0).
//! Default off: `mip_decompose`.

use crate::lp_data::highs_solve::highs_status_from_highs_model_status;
use crate::lp_data::presolve_support::get_kkt_failures;
use crate::prelude::*;

/// What one block's MIP returned.
struct BlockResult {
    status: HighsModelStatus,
    obj: Option<f64>,
    bound: f64,
    nodes: i64,
    iters: i64,
    col_value: Vec<f64>,
}

impl Highs {
    /// Components of the model: `(row_count, col_label, row_label)` with label = component
    /// index; columns without rows share one component.
    fn mip_components(&self) -> Option<(usize, Vec<usize>, Vec<usize>)> {
        let lp = &self.model_lp_;
        if lp.a_matrix_.format_ != MatrixFormat::kColwise {
            return None;
        }
        let (n, m) = (lp.num_col_ as usize, lp.num_row_ as usize);
        let mut par: Vec<usize> = (0..n + m).collect();
        fn find(par: &mut [usize], mut x: usize) -> usize {
            while par[x] != x {
                par[x] = par[par[x]];
                x = par[x];
            }
            x
        }
        let a = &lp.a_matrix_;
        for j in 0..n {
            for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                let (x, y) = (find(&mut par, j), find(&mut par, n + a.index_[p] as usize));
                if x != y {
                    par[x] = y;
                }
            }
        }
        let mut label = std::collections::HashMap::new();
        let mut lab = vec![0usize; n + m];
        // columns that appear in no row go to one shared component (label 0 is reserved)
        let mut has_free = false;
        for j in 0..n {
            if a.start_[j] == a.start_[j + 1] {
                has_free = true;
            }
        }
        let mut next = usize::from(has_free);
        for x in 0..n + m {
            let free_col = x < n && a.start_[x] == a.start_[x + 1];
            lab[x] = if free_col {
                0
            } else {
                let r = find(&mut par, x);
                *label.entry(r).or_insert_with(|| {
                    next += 1;
                    next - 1
                })
            };
        }
        (next >= 2).then(|| (next, lab[..n].to_vec(), lab[n..].to_vec()))
    }

    /// Solve each component as its own MIP and combine. `None` when the lever is off or the
    /// model does not split (the caller then runs the usual path).
    pub(crate) fn call_solve_mip_decomposed(&mut self) -> Option<HighsStatus> {
        if !self.options_.mip_decompose
            || self.solution_.value_valid
            || self.model_lp_.has_semi_variables()
        {
            return None;
        }
        let (k, cl, rl) = self.mip_components()?;
        let lp = &self.model_lp_;
        let (n, m) = (lp.num_col_ as usize, lp.num_row_ as usize);
        let mut cols: Vec<Vec<usize>> = vec![Vec::new(); k];
        let mut rows: Vec<Vec<usize>> = vec![Vec::new(); k];
        for j in 0..n {
            cols[cl[j]].push(j);
        }
        for i in 0..m {
            rows[rl[i]].push(i);
        }
        crate::highs_log_user!(
            &self.options_.log_options,
            HighsLogType::kInfo,
            "MIP decomposes into {} independent components\n",
            k
        );
        let mut order: Vec<usize> = (0..k).collect();
        // easiest first: the time of fast blocks flows to the big ones
        order.sort_by_key(|&cc| cols[cc].len() + rows[cc].len());
        let mut subs: Vec<HighsLp> = Vec::with_capacity(k);
        for &cc in &order {
            let a = &lp.a_matrix_;
            let mut rmap = vec![usize::MAX; m];
            for (i, &r) in rows[cc].iter().enumerate() {
                rmap[r] = i;
            }
            let mut sub = HighsLp::default();
            sub.num_col_ = cols[cc].len() as HighsInt;
            sub.num_row_ = rows[cc].len() as HighsInt;
            sub.sense_ = lp.sense_;
            sub.col_cost_ = cols[cc].iter().map(|&j| lp.col_cost_[j]).collect();
            sub.col_lower_ = cols[cc].iter().map(|&j| lp.col_lower_[j]).collect();
            sub.col_upper_ = cols[cc].iter().map(|&j| lp.col_upper_[j]).collect();
            sub.row_lower_ = rows[cc].iter().map(|&i| lp.row_lower_[i]).collect();
            sub.row_upper_ = rows[cc].iter().map(|&i| lp.row_upper_[i]).collect();
            if !lp.integrality_.is_empty() {
                sub.integrality_ = cols[cc].iter().map(|&j| lp.integrality_[j]).collect();
            }
            sub.a_matrix_.format_ = MatrixFormat::kColwise;
            sub.a_matrix_.num_col_ = sub.num_col_;
            sub.a_matrix_.num_row_ = sub.num_row_;
            sub.a_matrix_.start_ = vec![0];
            for &j in &cols[cc] {
                for p in a.start_[j] as usize..a.start_[j + 1] as usize {
                    sub.a_matrix_
                        .index_
                        .push(rmap[a.index_[p] as usize] as HighsInt);
                    sub.a_matrix_.value_.push(a.value_[p]);
                }
                sub.a_matrix_
                    .start_
                    .push(sub.a_matrix_.index_.len() as HighsInt);
            }
            subs.push(sub);
        }
        // Blocks run on `workers` threads (one per block up to `threads`). A block starting
        // when `left` seconds remain gets left * min(workers, blocks_left) / blocks_left, so
        // the time of finished blocks flows to the later ones. A block that hit its share
        // (kTimeLimit) is run again in a later pass with the time that is left, from its own
        // incumbent; its dual bound is the best over passes. For a sum of independent blocks
        // the objective cutoff of one block (global incumbent minus the others' bounds) is
        // never tighter than that block's own incumbent, so the start solution carries it.
        let workers = (self.options_.threads.max(1) as usize).min(k);
        let mut base_opts = self.options_.clone();
        base_opts.mip_decompose = false;
        if workers > 1 {
            base_opts.threads = 1;
            // interleaved child logs are unreadable; the per-block summary lines remain
            base_opts.output_flag = false;
            base_opts.log_options.output_flag = false;
        }
        let limit = self.options_.time_limit;
        let t_start = self.timer_.read();
        let clock = std::time::Instant::now();
        let stop = std::sync::atomic::AtomicBool::new(false);
        let best: std::sync::Mutex<Vec<Option<BlockResult>>> =
            std::sync::Mutex::new((0..k).map(|_| None).collect());
        let maximize = lp.sense_ == ObjSense::kMaximize;
        let mut pass = 0usize;
        loop {
            let pending: Vec<usize> = {
                let b = best.lock().unwrap_or_else(|e| e.into_inner());
                (0..k)
                    .filter(|&i| {
                        b[i].as_ref().is_none_or(|r| {
                            matches!(
                                r.status,
                                HighsModelStatus::kTimeLimit | HighsModelStatus::kInterrupt
                            )
                        })
                    })
                    .collect()
            };
            if pending.is_empty() || stop.load(std::sync::atomic::Ordering::Relaxed) {
                break;
            }
            if pass > 0
                && (!limit.is_finite() || limit - t_start - clock.elapsed().as_secs_f64() < 1.0)
            {
                break;
            }
            let np = pending.len();
            let next = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
            std::thread::scope(|scope| {
                for _ in 0..workers.min(np) {
                    scope.spawn(|| {
                        loop {
                            let slot = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if slot >= np || stop.load(std::sync::atomic::Ordering::Relaxed) {
                                return;
                            }
                            let idx = pending[slot];
                            let mut opts = base_opts.clone();
                            // The hard limit is all the time that is left; the share only
                            // binds (through the interrupt callback, checked in the tree)
                            // while blocks are still waiting to start.
                            let mut share = f64::INFINITY;
                            if limit.is_finite() {
                                let left =
                                    (limit - t_start - clock.elapsed().as_secs_f64()).max(0.0);
                                let left_blocks = (np - slot) as f64;
                                share = left * left_blocks.min(workers as f64) / left_blocks;
                                opts.time_limit = if workers > 1 { left } else { share };
                            }
                            let t_block = std::time::Instant::now();
                            let start_col = best.lock().unwrap_or_else(|e| e.into_inner())[idx]
                                .as_ref()
                                .filter(|r| r.obj.is_some())
                                .map(|r| r.col_value.clone());
                            let mut child = Highs::default();
                            child.pass_options(&opts);
                            child.pass_model(subs[idx].clone());
                            if share.is_finite() && workers > 1 {
                                let (next, t0) = (std::sync::Arc::clone(&next), t_block);
                                child.set_callback(HighsCallbackFunctionType::new(
                                    move |kind: HighsInt,
                                          _msg: &str,
                                          _out: &HighsCallbackOutput,
                                          cb_in: &mut HighsCallbackInput| {
                                        if kind == kCallbackMipInterrupt
                                            && next.load(std::sync::atomic::Ordering::Relaxed) < np
                                            && t0.elapsed().as_secs_f64() > share
                                        {
                                            cb_in.user_interrupt = true;
                                        }
                                    },
                                ));
                                child.start_callback(kCallbackMipInterrupt);
                            }
                            if let Some(cv) = start_col {
                                let sol = HighsSolution {
                                    col_value: cv,
                                    ..HighsSolution::default()
                                };
                                child.set_solution(&sol);
                            }
                            child.run();
                            let st = child.get_model_status();
                            let info = child.get_info();
                            let has_sol = child.solution_.value_valid
                                && info.objective_function_value.is_finite();
                            let mut res = BlockResult {
                                status: st,
                                obj: has_sol.then_some(info.objective_function_value),
                                bound: info.mip_dual_bound,
                                nodes: info.mip_node_count,
                                iters: i64::from(info.simplex_iteration_count.max(0)),
                                col_value: if has_sol {
                                    child.solution_.col_value.clone()
                                } else {
                                    Vec::new()
                                },
                            };
                            crate::highs_log_user!(
                                &self.options_.log_options,
                                HighsLogType::kInfo,
                                "decompose pass {} block {} ({}x{}): {:?} obj {} bound {} \
                                 {:.1}s\n",
                                pass + 1,
                                idx,
                                subs[idx].num_row_,
                                subs[idx].num_col_,
                                st,
                                res.obj.unwrap_or(f64::NAN),
                                res.bound,
                                t_block.elapsed().as_secs_f64()
                            );
                            if matches!(
                                st,
                                HighsModelStatus::kInfeasible | HighsModelStatus::kUnbounded
                            ) {
                                stop.store(true, std::sync::atomic::Ordering::Relaxed);
                            }
                            let mut b = best.lock().unwrap_or_else(|e| e.into_inner());
                            if let Some(old) = b[idx].as_ref() {
                                // keep the better incumbent and the better bound of both passes
                                res.nodes += old.nodes;
                                res.iters += old.iters;
                                res.bound = if maximize {
                                    res.bound.min(old.bound)
                                } else {
                                    res.bound.max(old.bound)
                                };
                                if let Some(oo) = old.obj
                                    && res
                                        .obj
                                        .is_none_or(|ro| if maximize { oo > ro } else { oo < ro })
                                {
                                    res.obj = Some(oo);
                                    res.col_value = old.col_value.clone();
                                }
                            }
                            b[idx] = Some(res);
                        }
                    });
                }
            });
            pass += 1;
        }
        let results = best.into_inner().unwrap_or_else(|e| e.into_inner());
        let mut col_value = vec![0.0; n];
        let (mut obj, mut bound) = (lp.offset_, lp.offset_);
        let (mut nodes, mut iters) = (0i64, 0i64);
        let (mut all_sol, mut all_opt) = (true, true);
        let mut status = HighsModelStatus::kOptimal;
        for (idx, res) in results.iter().enumerate() {
            let Some(r) = res else {
                // not started (a sibling was infeasible or unbounded)
                all_opt = false;
                all_sol = false;
                bound = f64::NEG_INFINITY;
                continue;
            };
            nodes += r.nodes;
            iters += r.iters;
            if r.status != HighsModelStatus::kOptimal {
                all_opt = false;
                let worse = matches!(
                    r.status,
                    HighsModelStatus::kInfeasible
                        | HighsModelStatus::kUnbounded
                        | HighsModelStatus::kUnboundedOrInfeasible
                );
                if worse || status == HighsModelStatus::kOptimal {
                    status = r.status;
                }
            }
            bound += r.bound;
            match r.obj {
                Some(o) => {
                    obj += o;
                    for (i, &j) in cols[order[idx]].iter().enumerate() {
                        col_value[j] = r.col_value[i];
                    }
                }
                None => all_sol = false,
            }
        }
        self.model_status_ = if all_opt {
            HighsModelStatus::kOptimal
        } else {
            status
        };
        let mut return_status = highs_status_from_highs_model_status(self.model_status_);
        if all_sol {
            self.solution_.col_value = col_value;
            self.model_lp_.a_matrix_.product_quad(
                &mut self.solution_.row_value,
                &self.solution_.col_value,
                kDebugReportOff,
            );
            self.solution_.value_valid = true;
        }
        self.info_.objective_function_value = if all_sol { obj } else { kHighsInf };
        let ptol = self.options_.primal_feasibility_tolerance;
        self.options_.primal_feasibility_tolerance = self.options_.mip_feasibility_tolerance;
        get_kkt_failures(
            &self.options_,
            false,
            &self.model_lp_,
            &self.model_lp_.col_cost_,
            &self.solution_,
            &mut self.info_,
            false,
        );
        self.options_.primal_feasibility_tolerance = ptol;
        self.info_.mip_node_count = nodes;
        self.info_.mip_dual_bound = bound;
        self.info_.mip_gap = if all_sol {
            (obj - bound).abs() / obj.abs().max(1.0)
        } else {
            kHighsInf
        };
        self.info_.simplex_iteration_count = iters.min(i64::from(kHighsIInf)) as HighsInt;
        self.info_.valid = true;
        if self.model_status_ == HighsModelStatus::kOptimal {
            return_status = self.check_optimality("MIP");
        }
        Some(return_status)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two independent blocks plus a column in no row: the combined optimum equals the whole.
    fn model(decompose: bool) -> Highs {
        let mut h = Highs::default();
        h.set_option_value_bool("output_flag", false);
        h.set_option_value_bool("mip_decompose", decompose);
        // x y | z w | v
        for (cost, lo, hi) in [
            (1.0, 0.0, 1.0),
            (2.0, 0.0, 1.0),
            (3.0, 0.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 2.0, 5.0),
        ] {
            h.add_col(cost, lo, hi, 0, &[], &[]);
        }
        for j in 0..5 {
            h.change_col_integrality(j, HighsVarType::kInteger);
        }
        h.add_row(1.0, kHighsInf, 2, &[0, 1], &[1.0, 1.0]);
        h.add_row(1.0, kHighsInf, 2, &[2, 3], &[1.0, 1.0]);
        h
    }

    #[test]
    fn components_solved_separately_match_the_whole() {
        let mut whole = model(false);
        whole.run();
        let mut split = model(true);
        split.model_lp_.a_matrix_.ensure_colwise();
        assert_eq!(split.mip_components().map(|t| t.0), Some(3));
        split.run();
        assert_eq!(split.get_model_status(), HighsModelStatus::kOptimal);
        assert_eq!(whole.get_model_status(), HighsModelStatus::kOptimal);
        let (a, b) = (
            whole.get_info().objective_function_value,
            split.get_info().objective_function_value,
        );
        assert!((a - 4.0).abs() < 1e-9 && (a - b).abs() < 1e-9, "{a} {b}");
        assert_eq!(split.solution_.col_value.len(), 5);
    }
}
