//! `rusolve solve --ekk` (or `RUSOLVE_EKK=1`): the LP solved by
//! `rusolve-ekk`, the faithful port of the HiGHS serial simplex, through the
//! ported `Highs::passModel` + `Highs::run` driver, with the ported LP
//! presolve by default (`RUSOLVE_EKK_PRESOLVE=off` disables it).
//!
//! The model is handed over as the `HighsLp` HiGHS's own free-format MPS
//! reader would build (MIT, `reference/HiGHS/highs/io/HMpsFF.cpp`), so that
//! the port reproduces HiGHS pivot for pivot:
//! - matrix columns keep their nonzeros in FILE order (`fillMatrix`,
//!   :136-170), zero values are skipped and the first of duplicate entries
//!   in a column wins (`parseColumns`, :911-956); the same rules for costs;
//! - row bounds from the row type and RHS (`parseRows` :630-642, `addRhs`
//!   :1065-1078), then RANGES (`parseRanges` `addRhs`, :1554-1565):
//!   L: lower = upper - |R|; G: upper = lower + |R|; E: R > 0 gives
//!   upper = lower + R, R < 0 gives lower = upper - |R|;
//! - integrality dropped (the LP relaxation, as `Highs::clearIntegrality`).
//!
//! Known gap: rusolve-io stores RANGES as a dense vector defaulting to 0, so
//! an explicit `R = 0` on an L/G row (which HiGHS turns into an equality)
//! reads as "no range". No Netlib / LP-47 file has one.

use rusolve_ekk::lp_data::highs_run::Highs;
use rusolve_ekk::prelude::*;
use rusolve_io::{Model, MpsFileOrder, ObjSense as IoSense, RowSense};
use std::process::ExitCode;

/// The `HighsLp` HiGHS reads from the same MPS file (see the module docs).
pub fn highs_lp_from_model(m: &Model, order: &MpsFileOrder) -> HighsLp {
    let n = m.num_variables;
    let nr = m.num_constraints;
    // Costs: first nonzero entry per column wins.
    let mut col_cost = vec![0.0; n];
    let mut has_cost = vec![false; n];
    for &(j, v) in &order.cost_entries {
        if v != 0.0 && !has_cost[j] {
            col_cost[j] = v;
            has_cost[j] = true;
        }
    }
    // Matrix: stable bucket by column (file order within a column), zeros
    // skipped, first of duplicate (row, col) entries kept.
    let mut count = vec![0usize; n + 1];
    for &(_, j, v) in &order.entries {
        if v != 0.0 {
            count[j + 1] += 1;
        }
    }
    for j in 0..n {
        count[j + 1] += count[j];
    }
    let mut next = count.clone();
    let total = count[n];
    let mut index = vec![0 as HighsInt; total];
    let mut value = vec![0.0; total];
    for &(i, j, v) in &order.entries {
        if v != 0.0 {
            index[next[j]] = i as HighsInt;
            value[next[j]] = v;
            next[j] += 1;
        }
    }
    let mut start: Vec<HighsInt> = Vec::with_capacity(n + 1);
    let mut kept_index = Vec::with_capacity(total);
    let mut kept_value = Vec::with_capacity(total);
    let mut seen = vec![usize::MAX; nr];
    start.push(0);
    for j in 0..n {
        for k in count[j]..count[j + 1] {
            let i = index[k] as usize;
            if seen[i] != j {
                seen[i] = j;
                kept_index.push(index[k]);
                kept_value.push(value[k]);
            }
        }
        start.push(kept_index.len() as HighsInt);
    }
    // Row bounds.
    let inf = kHighsInf;
    let mut row_lower = vec![0.0; nr];
    let mut row_upper = vec![0.0; nr];
    for i in 0..nr {
        let b = m.rhs[i];
        let (mut lo, mut up) = match m.senses[i] {
            RowSense::GreaterEqual => (b, inf),
            RowSense::Equal => (b, b),
            RowSense::LessEqual => (-inf, b),
        };
        let r = m.ranges[i];
        if r != 0.0 {
            match m.senses[i] {
                RowSense::LessEqual => lo = up - r.abs(),
                RowSense::GreaterEqual => up = lo + r.abs(),
                RowSense::Equal if r < 0.0 => lo = up - r.abs(),
                RowSense::Equal => up = lo + r.abs(),
            }
        }
        row_lower[i] = lo;
        row_upper[i] = up;
    }
    let mut lp = HighsLp {
        num_col_: n as HighsInt,
        num_row_: nr as HighsInt,
        col_cost_: col_cost,
        col_lower_: m.bounds.iter().map(|b| b.lower).collect(),
        col_upper_: m.bounds.iter().map(|b| b.upper).collect(),
        row_lower_: row_lower,
        row_upper_: row_upper,
        sense_: match m.obj_sense {
            IoSense::Minimize => ObjSense::kMinimize,
            IoSense::Maximize => ObjSense::kMaximize,
        },
        offset_: m.obj_offset,
        ..HighsLp::default()
    };
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = n as HighsInt;
    lp.a_matrix_.num_row_ = nr as HighsInt;
    lp.a_matrix_.start_ = start;
    lp.a_matrix_.index_ = kept_index;
    lp.a_matrix_.value_ = kept_value;
    lp
}

/// The names HiGHS's MPS reader hands over with the model (see the
/// comment on each field), shared by every path that loads a model for the
/// port: `mip_highs::read_mip` and the LP `solve --ekk`.
pub(crate) fn fill_names(lp: &mut HighsLp, model: &Model) {
    lp.col_names_ = model.col_names.clone();
    lp.row_names_ = model.row_names.clone();
    // `lp.objective_name_` (HMpsFF.cpp:88, default "Objective" from :600; the
    // uniqueness repair of `findModelObjectiveName` at :128 is not ported and
    // cannot trigger: the reader has one objective row and no name clash with a
    // column is checked here).
    lp.objective_name_ = model
        .obj_name
        .clone()
        .unwrap_or_else(|| "Objective".to_string());
}

/// The `HighsHessian` HiGHS's MPS reader hands over with the model
/// (HMpsFF.cpp:119-126): the square column-wise Hessian from
/// `Model::fill_hessian` (`q_dim`, `q_format`, `q_start`, `q_index`,
/// `q_value`), cleared when it has no starts (an LP: `dim_` 0).
pub fn highs_hessian_from_model(m: &Model) -> HighsHessian {
    let h = m.fill_hessian();
    let mut hessian = HighsHessian {
        dim_: h.dim,
        format_: HessianFormat::kSquare,
        start_: h.start,
        index_: h.index,
        value_: h.value,
    };
    // hessian must have at least start_[0]=0 for the fictitious column 0
    if hessian.start_.is_empty() {
        hessian.clear();
    }
    hessian
}

/// The CLI status name: the old engine's `SimplexStatus` spelling where one
/// exists, so the harness parsers stay unchanged.
fn status_name(s: HighsModelStatus) -> String {
    match s {
        HighsModelStatus::kOptimal => "Optimal".into(),
        HighsModelStatus::kInfeasible => "Infeasible".into(),
        HighsModelStatus::kUnbounded => "Unbounded".into(),
        HighsModelStatus::kIterationLimit | HighsModelStatus::kTimeLimit => "IterationLimit".into(),
        other => format!("{other:?}"),
    }
}

/// `--format json` / `--solution` for [`solve_ekk_cli`].
#[derive(Clone, Copy)]
pub struct JsonOut {
    pub json: bool,
    pub include_solution: bool,
}

/// `rusolve solve --ekk`: read, solve with the HiGHS-port driver, print the
/// simplex output lines (iterations reported as phase 2: HEkk keeps one
/// count). Exit codes are the exit-code table (`mip_highs.rs`): 0 ran, 3 input, 4 internal.
pub fn solve_ekk_cli(
    path: &str,
    csv: bool,
    duals: bool,
    time_limit: Option<f64>,
    iteration_limit: Option<usize>,
    out: JsonOut,
    solution_file: Option<&str>,
) -> ExitCode {
    let t_start = std::time::Instant::now();
    // Every rusolve_io format; MPS keeps its file order, the others the CSC order.
    let read = rusolve_io::read_model_file_order(path);
    let (model, order) = match read {
        Ok(x) => x,
        Err(e) => return crate::diag::fail_read(out.json, path, &e.to_string()),
    };
    if let Err(e) = model.gc.reject("rusolve solve") {
        eprintln!("{e}");
        return ExitCode::from(3);
    }
    // Convex QP: a model with a Hessian is checked before any engine starts.
    // A model without one skips this block and the Hessian branch below entirely.
    let qp_assessment = if model.hessian.is_empty() {
        None
    } else {
        match crate::qp::assess(&model) {
            Ok(a) => Some(a),
            Err(r) => {
                eprintln!("rusolve solve: {}", r.message);
                if out.json && !crate::diag::json_deferred() {
                    println!("{}", r.json);
                }
                return ExitCode::from(3);
            }
        }
    };
    let mut lp = highs_lp_from_model(&model, &order);
    // HiGHS always has the names after reading a model
    fill_names(&mut lp, &model);
    // The cheap lint rules, on a copy, before the solve (never touches the search)
    let warn = !crate::adopt::no_warnings(&std::env::args().collect::<Vec<_>>());
    let mut diags = if warn {
        // `solve` relaxes integrality, but the big-M rule is about the model as written
        if model.bounds.iter().any(|b| b.integer) {
            let mut as_written = lp.clone();
            as_written.integrality_ = model
                .bounds
                .iter()
                .map(|b| {
                    if b.integer {
                        HighsVarType::kInteger
                    } else {
                        HighsVarType::kContinuous
                    }
                })
                .collect();
            crate::lint::warnings(&as_written)
        } else {
            crate::lint::warnings(&lp)
        }
    } else {
        Vec::new()
    };
    crate::adopt::print_warnings(&diags);
    let mut options = HighsOptions {
        // RUSOLVE_EKK_LOG=1 keeps the HiGHS log on the console.
        output_flag: std::env::var("RUSOLVE_EKK_LOG").is_ok_and(|v| v != "0"),
        // HiGHS default ("choose") runs the ported LP presolve;
        // RUSOLVE_EKK_PRESOLVE=off keeps the bit-exact presolve-off path.
        presolve: std::env::var("RUSOLVE_EKK_PRESOLVE")
            .unwrap_or_else(|_| HighsOptions::default().presolve),
        ..crate::base_options()
    };
    if let Some(t) = time_limit {
        // `--time-limit` bounds the whole run, as in the legacy path (its deadline is
        // `t_start + t`): the model read (5 s on dlr1's 920 MB) is charged to it.
        options.time_limit = (t - t_start.elapsed().as_secs_f64()).max(0.0);
    }
    if let Some(k) = iteration_limit {
        options.simplex_iteration_limit = k.min(HighsInt::MAX as usize) as HighsInt;
    }
    // RUSOLVE_EKK_HIPO=1: the HiPO interior point code is reachable for
    // `solver=ipm` (the build switch of HiGHS 1.15+ with HiPO compiled in).
    if std::env::var("RUSOLVE_EKK_HIPO").is_ok_and(|v| v != "0") {
        options.hipo_available = true;
    }
    options.sync_log_options();
    let mut highs = Highs::with_options(options);
    // RUSOLVE_EKK_OPTS="key=value,...": HiGHS options by name (the LP A/B of
    // the default-off mip_* options that reach LP presolve).
    if let Ok(list) = std::env::var("RUSOLVE_EKK_OPTS") {
        for kv in list.split(',').filter(|s| !s.is_empty()) {
            let ok = kv
                .split_once('=')
                .is_some_and(|(k, v)| highs.set_option_value_str(k, v) != HighsStatus::kError);
            if !ok {
                eprintln!("RUSOLVE_EKK_OPTS: unknown option or invalid value: {kv}");
                return ExitCode::from(2);
            }
        }
    }
    if highs.pass_model(lp) == HighsStatus::kError {
        eprintln!("ekk: passModel failed");
        return ExitCode::from(4);
    }
    // HiGHS app: readModel (passModel(HighsModel), Highs.cpp:443-523) then
    // run. The model's Hessian is passed before run; an LP (dim 0) has none,
    // and passModel does no Hessian work for it (Highs.cpp:500, 512).
    let hessian = match &qp_assessment {
        Some(a) => crate::qp::solver_hessian(&model, a),
        None => highs_hessian_from_model(&model),
    };
    if hessian.dim_ > 0 && highs.pass_hessian(hessian) == HighsStatus::kError {
        eprintln!("ekk: passHessian failed");
        return ExitCode::from(4);
    }
    // Convex QP engine choice: with no explicit `solver`, a QP goes to HiPO first and
    // falls back to QUASS when HiPO does not end Optimal. Measured on a Maros-Meszaros slice:
    // HiPO solved 9 of 10 Maros-Meszaros instances in at most 0.2 s (AUG3DCQP 0.08 s against
    // 21 s for QUASS) and failed LISWET1 (QUASS 1.2 s). Naming `solver` keeps that choice.
    // Only QPs reach this; an LP or MIP never sets `qp_route`.
    let mut qp_route = qp_assessment.as_ref().map(|_| {
        let explicit = highs.get_options().solver != kHighsChooseString;
        if !explicit {
            highs.options_.solver = kIpmString.to_string();
            highs.options_.hipo_available = true;
        }
        crate::qp::QpRoute::new(explicit)
    });
    let mut run_status = highs.run();
    if let Some(route) = qp_route.as_mut()
        && !route.explicit
    {
        let first = highs.get_model_status();
        let spent_out = matches!(
            first,
            HighsModelStatus::kTimeLimit | HighsModelStatus::kInterrupt
        );
        if first != HighsModelStatus::kOptimal && !spent_out {
            route.fallback_from = Some(status_name(first));
            highs.options_.solver = kSimplexString.to_string();
            if let Some(t) = time_limit {
                highs.options_.time_limit = (t - t_start.elapsed().as_secs_f64()).max(0.0);
            }
            run_status = highs.run();
        }
    }
    let status = highs.get_model_status();
    let info = highs.get_info();
    let sol = highs.get_solution();
    let name = status_name(status);
    let obj = info.objective_function_value;
    let iters = info.simplex_iteration_count.max(0);
    // Infeasible set, ray, violated row (after the solve, bounded budget)
    let time_s = t_start.elapsed().as_secs_f64();
    let post = crate::adopt::post_solve(
        path,
        &highs.model_lp_,
        status,
        sol,
        warn,
        crate::adopt::budget_for(time_s, time_limit.unwrap_or(f64::INFINITY)),
    );
    crate::adopt::print_warnings(&post);
    diags.extend(post);
    if out.json {
        // The envelope the port MIP path prints (`envelope.rs`), from the same `Highs`.
        let solved = crate::envelope::Solved {
            path,
            lp: &highs.model_lp_,
            status,
            run_status,
            info,
            solution: sol,
            time_s,
            seed: highs.get_options().random_seed,
            threads: highs.get_options().threads,
            time_limit: highs.get_options().time_limit,
            node_limit: None,
            mip_max_nodes: highs.get_options().mip_max_nodes,
            options: &[],
            include_solution: out.include_solution,
            // an LP solve: no MIP statistics, no hook component
            mip_stats: None,
            wasm: None,
        };
        let text = crate::adopt::with_diagnostics(&crate::envelope::envelope_json(&solved), &diags);
        match &qp_assessment {
            None => println!("{text}"),
            Some(a) => {
                let route = qp_route.as_ref().expect("a QP has a route");
                println!("{}", crate::qp::with_qp(&text, a, route));
            }
        }
    } else if csv {
        println!("status,objective,phase1_iters,phase2_iters,x");
        let xs: Vec<String> = sol.col_value.iter().map(|v| format!("{v}")).collect();
        println!("{name},{obj},0,{iters},[{}]", xs.join(";"));
    } else {
        println!("status:  {name}");
        println!("objective: {obj}");
        println!("iterations: 0 (phase1) + {iters} (phase2) [ekk]");
        if let (Some(a), Some(route)) = (&qp_assessment, &qp_route) {
            println!("{}", crate::qp::text_line(a, route));
        }
        let first5: Vec<String> = sol
            .col_value
            .iter()
            .take(5)
            .map(|v| format!("{v:.6}"))
            .collect();
        println!("x[:5]: [{}]", first5.join(", "));
    }
    // `--solution-file PATH`: the HiGHS raw solution file with the duals the engine produced,
    // which `rusolve check` reads for the QP KKT residuals.
    if let Some(p) = solution_file {
        let lo = HighsLogOptions {
            log_to_console: false,
            ..Default::default()
        };
        if let Err(e) = crate::envelope::write_solution_file(
            p,
            lo,
            &crate::envelope::Solved {
                path,
                lp: &highs.model_lp_,
                status,
                run_status,
                info,
                solution: sol,
                time_s: 0.0,
                seed: 0,
                threads: 1,
                time_limit: f64::INFINITY,
                node_limit: None,
                mip_max_nodes: 0,
                options: &[],
                include_solution: false,
                mip_stats: None,
                wasm: None,
            },
        ) {
            eprintln!("rusolve solve: cannot write {p}: {e}");
            return ExitCode::from(4);
        }
    }
    // stderr is for warnings and errors: a QP that ended Optimal prints nothing there (its class,
    // engine and any fallback are in the text line / the `qp` block of the envelope).
    if qp_assessment.is_none()
        || run_status != HighsStatus::kOk
        || status != HighsModelStatus::kOptimal
    {
        eprintln!("rusolve ekk: run_status {run_status:?} model_status {status:?}");
    }
    if duals {
        if sol.dual_valid {
            println!("row_duals: {:?}", sol.row_dual);
            println!("reduced_costs: {:?}", sol.col_dual);
        } else {
            println!("duals: unavailable (non-optimal status)");
        }
    }
    // The solve ran; the status is in the output, not in the exit code.
    if run_status == HighsStatus::kError {
        ExitCode::from(4)
    } else {
        ExitCode::SUCCESS
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MPS: &str = "\
NAME T
ROWS
 N obj
 L r1
 G r2
 E r3
 E r4
COLUMNS
 x obj 1 r3 1
 x r1 0.5 r1 9
 x r2 0
 y obj 0 obj 5
 y r4 1 r2 1
RHS
 rhs r1 4 r2 1
 rhs r3 6 r4 2
RANGES
 rng r1 1.5 r2 -2
 rng r3 3 r4 -1
ENDATA
";

    #[test]
    fn highs_mps_semantics() {
        let (m, o) = rusolve_io::read_mps_file_order(MPS).unwrap();
        let lp = highs_lp_from_model(&m, &o);
        // file order kept, zero skipped, first duplicate (x, r1) = 0.5 wins
        assert_eq!(lp.a_matrix_.start_, vec![0, 2, 4]);
        assert_eq!(lp.a_matrix_.index_, vec![2, 0, 3, 1]);
        assert_eq!(lp.a_matrix_.value_, vec![1.0, 0.5, 1.0, 1.0]);
        // zero cost skipped, then the first nonzero wins
        assert_eq!(lp.col_cost_, vec![1.0, 5.0]);
        let inf = kHighsInf;
        assert_eq!(lp.row_lower_, vec![2.5, 1.0, 6.0, 1.0]);
        assert_eq!(lp.row_upper_, vec![4.0, 3.0, 9.0, 2.0]);
        assert_eq!(lp.col_lower_, vec![0.0, 0.0]);
        assert_eq!(lp.col_upper_, vec![inf, inf]);
    }

    #[test]
    fn ekk_cli_solves_small_lp() {
        let (m, o) = rusolve_io::read_mps_file_order(MPS).unwrap();
        let mut options = HighsOptions {
            output_flag: false,
            presolve: kHighsOffString.to_string(),
            ..HighsOptions::default()
        };
        options.sync_log_options();
        let mut highs = Highs::with_options(options);
        assert_eq!(
            highs.pass_model(highs_lp_from_model(&m, &o)),
            HighsStatus::kOk
        );
        highs.run();
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        // x in [6, 8] (r3, r1), y in [1, 2] (r2, r4): min x + 5y = 11
        assert_eq!(highs.get_info().objective_function_value, 11.0);
    }

    const QP_MPS: &str = "\
NAME QP
ROWS
 N obj
 L r1
COLUMNS
 x obj -2 r1 1
 y obj -4 r1 1
RHS
 rhs r1 10
QUADOBJ
 x x 2
 y y 4
ENDATA
";

    #[test]
    fn qp_hessian_from_mps_is_square_csc() {
        let (m, _) = rusolve_io::read_mps_file_order(QP_MPS).unwrap();
        let h = highs_hessian_from_model(&m);
        assert_eq!(h.dim_, 2);
        assert_eq!(h.format_, HessianFormat::kSquare);
        assert_eq!(h.start_, vec![0, 1, 2]);
        assert_eq!(h.index_, vec![0, 1]);
        assert_eq!(h.value_, vec![2.0, 4.0]);
    }

    #[test]
    fn qp_lp_has_no_hessian() {
        let (m, _) = rusolve_io::read_mps_file_order(MPS).unwrap();
        // HMpsFF.cpp:126: no starts, so hessian.clear()
        let h = highs_hessian_from_model(&m);
        assert_eq!(h.dim_, 0);
        assert_eq!(h.format_, HessianFormat::kTriangular);
        assert_eq!(h.start_, vec![0]);
        assert!(h.index_.is_empty() && h.value_.is_empty());
    }

    #[test]
    fn qp_mps_qp_solves_with_passed_hessian() {
        // min 0.5 (2x^2 + 4y^2) - 2x - 4y, x + y <= 10, x, y >= 0: x = 1, y = 1
        let (m, o) = rusolve_io::read_mps_file_order(QP_MPS).unwrap();
        let mut options = HighsOptions {
            output_flag: false,
            presolve: kHighsOffString.to_string(),
            ..HighsOptions::default()
        };
        options.sync_log_options();
        let mut highs = Highs::with_options(options);
        assert_eq!(
            highs.pass_model(highs_lp_from_model(&m, &o)),
            HighsStatus::kOk
        );
        assert_eq!(
            highs.pass_hessian(highs_hessian_from_model(&m)),
            HighsStatus::kOk
        );
        assert_eq!(highs.run(), HighsStatus::kOk);
        assert_eq!(highs.get_model_status(), HighsModelStatus::kOptimal);
        let x = &highs.get_solution().col_value;
        // QUASS adds qp_regularization_value (1e-7) to the Hessian diagonal
        // (a_quass.cpp:149), so HiGHS returns 2/(2+1e-7) and 4/(4+1e-7).
        let reg = 1e-7;
        assert!((x[0] - 2.0 / (2.0 + reg)).abs() <= 1e-12);
        assert!((x[1] - 4.0 / (4.0 + reg)).abs() <= 1e-12);
        assert!((highs.get_info().objective_function_value + 3.0).abs() <= 1e-6);
    }
}
