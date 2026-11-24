//! Convex MIQP by LP-based outer approximation: the glue between the model files, the port's
//! MIP search and `rusolve_ekk::mip::oa` (which holds the cuts, the separator and the exactness
//! oracle, and cites the sources).
//!
//! Only `rusolve mip` / `rusolve solve` with `miqp_oa` on and a proved-convex Q get here
//! (`route.rs`); a plain LP or MIP never does. The flow, called from `mip_highs::run_mip_highs`:
//!
//! 1. [`prepare`]: read Q, refuse a Q that is not positive semidefinite (the cuts would be
//!    invalid), keep the user's LP, append the epigraph columns to the LP that goes to the port;
//! 2. [`run`]: solve the continuous QP relaxation (HiPO, QUASS on failure), add its tangent as
//!    rows, run the tree separator and the exactness loop;
//! 3. [`finish`]: recompute the incumbent's objective exactly with Q, cut the solution back to
//!    the user's columns and rows, and make the status honest (an unproved candidate is never
//!    `Optimal`).

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::mip::oa::{
    OaConfig, OaForm, OaOutcome, OaStats, OaStatsHandle, QuadForm, add_root_tangents, run_oa,
};
use rusolve_ekk::prelude::*;
use rusolve_io::Model;
use serde_json::{Value, json};
use std::sync::{Arc, Mutex};
use std::time::Instant;

/// What `prepare` keeps for the solve and the report.
pub struct OaPlan {
    /// The user's model, before the epigraph columns.
    pub user_lp: HighsLp,
    model: Model,
    qf: QuadForm,
    form: Arc<OaForm>,
    glo: Arc<Vec<f64>>,
    ghi: Arc<Vec<f64>>,
    pub tree: bool,
    pub disaggregate: bool,
    stats: OaStatsHandle,
    root: RootInfo,
}

#[derive(Default, Clone)]
struct RootInfo {
    status: String,
    objective: Option<f64>,
    tangents: usize,
}

/// Reads Q and appends the epigraph columns to `lp`. `Err` is the message of a typed refusal.
pub fn prepare(
    path: &str,
    lp: &mut HighsLp,
    disaggregate: bool,
    tree: bool,
) -> Result<OaPlan, String> {
    let (model, _) =
        rusolve_io::read_model_file_order(path).map_err(|e| format!("read {path}: {e}"))?;
    if let Err(why) = crate::qp::assess_miqp(&model) {
        return Err(format!("miqp_oa needs a convex quadratic objective: {why}"));
    }
    let n = lp.num_col_ as usize;
    let maximize = lp.sense_ == ObjSense::kMaximize;
    let qf = QuadForm {
        n,
        entries: model.hessian.clone(),
    };
    let form = OaForm::new(&qf, maximize, disaggregate);
    let user_lp = lp.clone();
    form.extend_lp(lp);
    Ok(OaPlan {
        glo: Arc::new(user_lp.col_lower_.clone()),
        ghi: Arc::new(user_lp.col_upper_.clone()),
        user_lp,
        model,
        qf,
        form: Arc::new(form),
        tree,
        disaggregate,
        stats: Arc::new(Mutex::new(OaStats::default())),
        root: RootInfo::default(),
    })
}

/// The continuous QP relaxation: HiPO first, QUASS when HiPO does not end optimal (the same
/// choice as the convex QP route, `ekk.rs`).
fn solve_root_qp(plan: &OaPlan, deadline: Option<Instant>) -> (HighsModelStatus, Vec<f64>, f64) {
    let mut rl = plan.user_lp.clone();
    rl.integrality_.clear();
    let mut options = crate::base_options();
    options.output_flag = false;
    options.threads = 1;
    let left = |d: Option<Instant>| {
        d.map_or(kHighsInf, |d| {
            d.saturating_duration_since(Instant::now()).as_secs_f64()
        })
    };
    options.time_limit = left(deadline);
    options.solver = kIpmString.to_string();
    options.hipo_available = true;
    options.sync_log_options();
    let mut h = Highs::with_options(options);
    h.pass_model(rl);
    h.pass_hessian(crate::ekk::highs_hessian_from_model(&plan.model));
    h.run();
    let mut status = h.get_model_status();
    if status != HighsModelStatus::kOptimal
        && !matches!(
            status,
            HighsModelStatus::kTimeLimit | HighsModelStatus::kInterrupt
        )
    {
        h.options_.solver = kSimplexString.to_string();
        h.options_.time_limit = left(deadline);
        h.run();
        status = h.get_model_status();
    }
    let obj = h.get_info().objective_function_value;
    (status, h.solution_.col_value.clone(), obj)
}

/// Root QP, root tangent, then the outer-approximation search on the loaded OA model.
pub fn run(plan: &mut OaPlan, highs: &mut Highs) -> OaOutcome {
    let tl = highs.options_.time_limit;
    let deadline =
        (tl.is_finite()).then(|| Instant::now() + std::time::Duration::from_secs_f64(tl));
    let (status, x, obj) = solve_root_qp(plan, deadline);
    plan.root.status = format!("{status:?}");
    if status == HighsModelStatus::kOptimal && x.len() == plan.form.n {
        plan.root.objective = Some(obj);
        plan.root.tangents =
            add_root_tangents(highs, &plan.form, &x, &plan.glo, &plan.ghi, &plan.stats);
    }
    if plan.tree {
        // hooks run in the presolved column space: the epigraph columns and the cut indices
        // are the model's own only with presolve off
        highs.set_option_value_str("presolve", "off");
    }
    run_oa(
        highs,
        Arc::clone(&plan.form),
        Arc::clone(&plan.glo),
        Arc::clone(&plan.ghi),
        Arc::clone(&plan.stats),
        OaConfig {
            tree: plan.tree,
            ..OaConfig::default()
        },
        deadline,
    )
}

/// After the search: exact objective, honest status, the solution cut to the user's size, and
/// the `miqp_oa` object of the envelope.
pub fn finish(plan: &OaPlan, highs: &mut Highs, out: &OaOutcome) -> Value {
    let n = plan.user_lp.num_col_ as usize;
    let m = plan.user_lp.num_row_ as usize;
    let nc = plan.form.num_cols();
    let tol = 2.0 * highs.options_.mip_feasibility_tolerance.max(1e-9);
    let mut status = highs.get_model_status();
    let mut max_viol = 0.0f64;
    let mut exact_obj = None;
    let has_solution = highs.solution_.value_valid && highs.solution_.col_value.len() >= nc;
    if has_solution {
        let z = &highs.solution_.col_value;
        let x = &z[..n];
        let obj = plan.user_lp.offset_
            + plan
                .user_lp
                .col_cost_
                .iter()
                .zip(x)
                .map(|(cost, v)| cost * v)
                .sum::<f64>()
            + plan.qf.value(x);
        exact_obj = Some(obj);
        // how far each epigraph value is below (above, for a maximisation) the quadratic
        for (k, b) in plan.form.blocks.iter().enumerate() {
            let sub = QuadForm {
                n,
                entries: block_entries(&plan.qf, &plan.form, k),
            };
            let f = sub.value(x);
            max_viol = max_viol.max(plan.form.sigma * (f - z[b.eta]));
        }
        let info = &mut highs.info_;
        info.objective_function_value = obj;
        let db = info.mip_dual_bound;
        info.mip_gap = if db.is_finite() {
            if obj != 0.0 {
                (obj - db).abs() / obj.abs()
            } else if (obj - db).abs() < 1e-9 {
                0.0
            } else {
                f64::INFINITY
            }
        } else {
            f64::INFINITY
        };
    }
    // An optimal OA model whose last candidate was never accepted (round limit) proved nothing.
    if status == HighsModelStatus::kOptimal && !out.lazy.accepted {
        status = HighsModelStatus::kIterationLimit;
    }
    // Without a tangent the first LP can be unbounded for a bounded MIQP: do not claim it.
    if matches!(
        status,
        HighsModelStatus::kUnbounded | HighsModelStatus::kUnboundedOrInfeasible
    ) && plan.root.status != "kUnbounded"
        && plan.root.status != "kUnboundedOrInfeasible"
    {
        status = HighsModelStatus::kUnknown;
    }
    highs.model_status_ = status;
    let s = &mut highs.solution_;
    s.col_value.truncate(n);
    s.row_value.truncate(m);
    s.col_dual.truncate(n);
    s.row_dual.truncate(m);
    let st = &out.stats;
    json!({
        "engine": "oa",
        "blocks": plan.form.blocks.len(),
        "disaggregated": plan.disaggregate,
        "tree_separation": plan.tree,
        "rounds": out.lazy.rounds,
        "accepted": out.lazy.accepted,
        "nodes_total": out.lazy.mip_nodes,
        "cuts": {
            "root_tangent": st.root_tangents,
            "tree": st.tree_cuts,
            "tree_calls": st.tree_calls,
            "lazy_new": st.lazy_cuts_new,
            "lazy_from_tree": st.lazy_cuts_from_pool,
            "dropped_unsafe": st.cuts_unsafe,
            "rows_added_to_model": out.lazy.rows_added,
        },
        "root_qp": {
            "status": plan.root.status,
            "objective": plan.root.objective,
        },
        "exact_objective": exact_obj,
        "max_epigraph_violation": has_solution.then_some(max_viol),
        "violation_tolerance": tol,
    })
}

/// Entries of `q` that fall in block `k` (for the exact value of one block).
fn block_entries(q: &QuadForm, form: &OaForm, k: usize) -> Vec<(usize, usize, f64)> {
    let cols = &form.blocks[k].cols;
    q.entries
        .iter()
        .copied()
        .filter(|&(i, _, _)| cols.binary_search(&i).is_ok())
        .collect()
}
