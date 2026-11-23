//! Golden test for the IPX linear algebra (iterate, KKT solvers,
//! conjugate residuals, normal matrices, preconditioners, power method,
//! basis, LU wrappers, crash basis, maxvolume, symbolic invert, starting
//! basis, model postsolve).
//!
//! Fixtures: `ipx/linalg/<lp>.<cfg>.jsonl`, written by `highs-oracle
//! ipx-linalg` (tools/highs-oracle/src/ipx.cpp): the real C++ IPM (
//! not replayed here) produces iterates, and the linear-algebra components are dumped on
//! them. Configurations: `default` (BasicLu, maxvolume heuristic),
//! `lu1` (`lu_kernel=1` ForrestTomlin, `update_heuristic=0` sequential
//! maxvolume), `dual` (`dualize=1`) and `slack` (`crash_basis=0`).
//! Committed: afiro.default, ipxmix.{default,lu1,dual}, ipxtall.default
//! (dualized by default), ipxdep.{default,dual,slack} (dependent rows,
//! duplicate columns); target-only: afiro, sc50a, kb2, adlittle, blend,
//! share2b, ipxmix, ipxtall, ipxdep in every configuration.
//!
//! Replay, bit for bit, starting each stage from the oracle state in the
//! record (iterates) or from the port's own state of the previous stage
//! (basis, iterate after StartingBasis / KKTSolverBasis::Factorize):
//! * `iterate`: Iterate construction (tolerances, bounds/costs measures) and
//!   the evaluation of a loaded state (objectives, residuals, mu, ...).
//! * `kkt_diag`: KKTSolverDiag Factorize (with and without iterate) + Solve
//!   (ConjugateResiduals on NormalMatrix with DiagonalPrecond).
//! * `starting_basis`: StartingBasis (CrashBasis, maxvolume, dependent
//!   rows/cols, LU factors and BasicLu memory state).
//! * `basis_ops`: SolveDense N/T/aliased, MinSingularValue (PowerMethod),
//!   SymbolicInvert, DensityInverse, RandomPermute.
//! * `kkt_basis`: KKTSolverBasis Factorize (DropPrimal/DropDual, maxvolume
//!   updates) + Solve (ConjugateResiduals on SplittedNormalMatrix).
//! * `basis_load`: Basis::Load of the statuses after a few main IPM
//!   iterations.
//! * `postsolve`: Model::PostsolveInteriorSolution, Iterate::Postprocess,
//!   EvaluatePostsolved, DropToComplementarity, ComputeBasicSolution,
//!   PostsolveBasicSolution, EvaluateBasicSolution, PostsolveBasis.

mod common;
mod port;

use std::rc::Rc;

use common::*;
use port::ipx::*;
use port::*;
use rusolve_ekk::ipm::ipx::basis::Basis;
use rusolve_ekk::ipm::ipx::control::Control;
use rusolve_ekk::ipm::ipx::iterate::Iterate;
use rusolve_ekk::ipm::ipx::kkt_solver::KKTSolver;
use rusolve_ekk::ipm::ipx::kkt_solver_basis::KKTSolverBasis;
use rusolve_ekk::ipm::ipx::kkt_solver_diag::KKTSolverDiag;
use rusolve_ekk::ipm::ipx::model::Model;
use rusolve_ekk::ipm::ipx::starting_basis::starting_basis;
use rusolve_ekk::ipm::ipx::symbolic_invert::random_permute_for_oracle;
use rusolve_ekk::ipm::ipx::{
    IPX_basic, IPX_nonbasic, IPX_nonbasic_lb, IPX_nonbasic_ub, IPX_superbasic, Info, Int, Vector,
};

struct Replay<'f> {
    fx: &'f Fixture,
    control: Rc<Control>,
    model: Rc<Model>,
    kkt_maxiter: Int,
}

impl Replay<'_> {
    fn iterate_rec(&self, tag: &str) -> Option<Rec<'_>> {
        self.fx.of_type("iterate").find(|r| r.str("tag") == tag)
    }

    fn load(&self, tag: &str) -> Iterate {
        let e = self.iterate_rec(tag).expect("iterate record");
        ipx_load_iterate(&self.control, &self.model, &e)
    }

    /// kktSolve of ipx.cpp: Factorize, then (no errflag) Solve of the
    /// record's random right-hand sides. Returns the info after the call.
    fn kkt(
        &self,
        chk: &mut Checker,
        e: &Rec,
        kkt: &mut dyn KKTSolver,
        mut it: Option<&mut Iterate>,
    ) -> Info {
        let mut info = ipx_info(&e.obj("info_in"));
        kkt.factorize(it.as_deref_mut(), &mut info);
        check_ipx_info(chk, "info_factorize", &e.obj("info_factorize"), &info);
        chk.int("basis_changes", e.int("basis_changes"), kkt.basis_changes());
        if let (Some(ei), Some(it)) = (e.opt_obj("iterate"), it.as_deref()) {
            check_ipx_iterate(chk, "iterate", &ei, it);
        }
        if let (Some(eb), Some(b)) = (e.opt_obj("basis"), kkt.basis()) {
            check_ipx_basis(chk, "basis", &eb, b);
        }
        if info.errflag != 0 || !e.has("a") {
            chk.bool("solved", e.has("a"), false);
            return info;
        }
        let (a, b) = (vector_of(e, "a"), vector_of(e, "b"));
        let (n, m) = (self.model.cols(), self.model.rows());
        let (mut x, mut y) = (Vector::new(n + m), Vector::new(m));
        kkt.solve(&a, &b, e.f64("tol"), &mut x, &mut y, &mut info);
        check_vec(chk, "solve", e, "x", &x);
        check_vec(chk, "solve", e, "y", &y);
        chk.int("iter", e.int("iter"), kkt.iter());
        check_ipx_info(chk, "info_solve", &e.obj("info_solve"), &info);
        info
    }
}

fn basis_ops(chk: &mut Checker, e: &Rec, basis: &Basis) {
    let m = basis.model().rows();
    let rhs = vector_of(e, "rhs");
    let (mut lhs_n, mut lhs_t) = (Vector::new(m), Vector::new(m));
    basis.solve_dense(&rhs, &mut lhs_n, b'N');
    basis.solve_dense(&rhs, &mut lhs_t, b'T');
    // SolveDense(alias, alias, 'N') as starting_basis.cc:76: the port
    // passes a copy of the right-hand side.
    let mut alias = rhs.clone();
    let tmp = alias.clone();
    basis.solve_dense(&tmp, &mut alias, b'N');
    check_vec(chk, "solve", e, "lhs_n", &lhs_n);
    check_vec(chk, "solve", e, "lhs_t", &lhs_t);
    check_vec(chk, "solve", e, "lhs_alias", &alias);
    chk.f64(
        "min_singular_value",
        e.f64("min_singular_value"),
        basis.min_singular_value(),
    );
    let (mut rowcounts, mut colcounts) = (vec![0; m as usize], vec![0; m as usize]);
    basis.symbolic_invert(Some(&mut rowcounts), Some(&mut colcounts));
    chk.ints("rowcounts", &e.ints("rowcounts"), &rowcounts);
    chk.ints("colcounts", &e.ints("colcounts"), &colcounts);
    chk.f64(
        "density_inverse",
        e.f64("density_inverse"),
        basis.density_inverse(),
    );
    chk.ints(
        "random_permute",
        &e.ints("random_permute"),
        &random_permute_for_oracle(basis.oracle_view().basis),
    );
}

fn postsolve(
    chk: &mut Checker,
    e: &Rec,
    model: &Model,
    it: &mut Iterate,
    basis: &Basis,
    info: &Info,
) {
    let (n, m) = (model.cols(), model.rows());
    let (nu, mu) = (e.f64s("x").len(), e.f64s("slack").len());
    it.postprocess();
    let mut pinfo = *info;
    it.evaluate_postsolved(&mut pinfo);
    check_ipx_info(chk, "info", &e.obj("info"), &pinfo);
    check_ipx_iterate(chk, "iterate", &e.obj("iterate"), it);
    let (mut xu, mut xlu, mut xuu, mut zlu, mut zuu) = (
        vec![0.0; nu],
        vec![0.0; nu],
        vec![0.0; nu],
        vec![0.0; nu],
        vec![0.0; nu],
    );
    let (mut su, mut yu) = (vec![0.0; mu], vec![0.0; mu]);
    model.postsolve_interior_solution(
        it.x(),
        it.xl(),
        it.xu(),
        it.y(),
        it.zl(),
        it.zu(),
        Some(&mut xu),
        Some(&mut xlu),
        Some(&mut xuu),
        Some(&mut su),
        Some(&mut yu),
        Some(&mut zlu),
        Some(&mut zuu),
    );
    for (k, a) in [
        ("x", &xu),
        ("xl", &xlu),
        ("xu", &xuu),
        ("slack", &su),
        ("y", &yu),
        ("zl", &zlu),
        ("zu", &zuu),
    ] {
        check_vec(chk, "interior", e, k, a);
    }
    let (mut xc, mut yc, mut zc) = (Vector::new(n + m), Vector::new(m), Vector::new(n + m));
    it.drop_to_complementarity(&mut xc, &mut yc, &mut zc);
    check_vec(chk, "crossover", e, "x_crossover", &xc);
    check_vec(chk, "crossover", e, "y_crossover", &yc);
    check_vec(chk, "crossover", e, "z_crossover", &zc);
    let basic_status: Vec<Int> = (0..n + m)
        .map(|j| {
            let x = xc[j as usize];
            if basis.is_basic(j) {
                IPX_basic
            } else if x == model.lb_at(j) {
                IPX_nonbasic_lb
            } else if x == model.ub_at(j) {
                IPX_nonbasic_ub
            } else if j < n {
                IPX_superbasic
            } else if model.lb_at(j).is_finite() {
                IPX_nonbasic_lb
            } else if model.ub_at(j).is_finite() {
                IPX_nonbasic_ub
            } else {
                IPX_nonbasic
            }
        })
        .collect();
    chk.ints("basic_status", &e.ints("basic_status"), &basic_status);
    let (mut xb, mut yb, mut zb) = (xc.clone(), yc.clone(), zc.clone());
    basis.compute_basic_solution(&mut xb, &mut yb, &mut zb);
    check_vec(chk, "basic", e, "x_basic", &xb);
    check_vec(chk, "basic", e, "y_basic", &yb);
    check_vec(chk, "basic", e, "z_basic", &zb);
    let (mut xbu, mut sbu, mut ybu, mut zbu) =
        (vec![0.0; nu], vec![0.0; mu], vec![0.0; mu], vec![0.0; nu]);
    model.postsolve_basic_solution(
        &xb,
        &yb,
        &zb,
        &basic_status,
        Some(&mut xbu),
        Some(&mut sbu),
        Some(&mut ybu),
        Some(&mut zbu),
    );
    check_vec(chk, "basic", e, "xb_user", &xbu);
    check_vec(chk, "basic", e, "slackb_user", &sbu);
    check_vec(chk, "basic", e, "yb_user", &ybu);
    check_vec(chk, "basic", e, "zb_user", &zbu);
    let mut binfo = Info::default();
    model.evaluate_basic_solution(&xb, &yb, &zb, &basic_status, &mut binfo);
    check_ipx_info(chk, "info_basic", &e.obj("info_basic"), &binfo);
    let (mut cbasis, mut vbasis) = (vec![0; mu], vec![0; nu]);
    model.postsolve_basis(&basic_status, Some(&mut cbasis), Some(&mut vbasis));
    chk.ints("cbasis", &e.ints("cbasis"), &cbasis);
    chk.ints("vbasis", &e.ints("vbasis"), &vbasis);
}

/// Counts of what the fixtures reached, so a regenerated fixture set that
/// silently stops exercising a component fails the test.
#[derive(Default)]
struct Reach {
    kkt_diag_solves: usize,
    kkt_basis_solves: usize,
    basis_changes: i64,
    dropped: i64,
    forrest_tomlin: usize,
    dependent: i64,
    dualized: usize,
    postsolves: usize,
}

fn replay(fails: &mut Failures, fx: &Fixture, reach: &mut Reach) {
    let control = Rc::new(ipx_control(fx));
    let input = fx.first("input");
    let em = fx.first("model");
    let mut model = Model::new();
    let err = ipx_load_model(&control, &input, &mut model);
    if err != 0 || em.int("errflag") != 0 {
        fails.run(&em.ctx, || {
            let mut chk = Checker::new(em.ctx.clone());
            chk.int("errflag", em.int("errflag"), err);
            Some(chk)
        });
        return;
    }
    reach.dualized += model.dualized() as usize;
    let m = model.rows();
    let rp = Replay {
        fx,
        control: control.clone(),
        model: Rc::new(model),
        kkt_maxiter: 500.min(10 + m / 20),
    };
    let model = rp.model.clone();

    // Port state carried from stage to stage.
    let mut it: Option<Iterate> = None;
    let mut basis: Option<Basis> = None;
    let mut basis2: Option<Basis> = None;
    let mut info = Info::default();
    for e in fx.recs().skip(3) {
        let tag = if e.has("tag") { e.str("tag") } else { "" };
        let ctx = format!("{}[{tag}]", e.ctx);
        fails.run(&ctx.clone(), || {
            let mut chk = Checker::new(ctx);
            match e.t() {
                "iterate" => {
                    let tag = e.str("tag");
                    let a = if tag == "constructed" {
                        ipx_new_iterate(&control, &model)
                    } else {
                        ipx_load_iterate(&control, &model, &e)
                    };
                    check_ipx_iterate(&mut chk, "iterate", &e, &a);
                }
                "ipm_starting_point" | "ipm_initial" | "ipm_main" => {}
                "kkt_diag" => {
                    let mut kkt = KKTSolverDiag::new(control.clone(), model.clone());
                    kkt.set_maxiter(rp.kkt_maxiter);
                    let tag = e.str("tag");
                    let mut a = (tag != "starting_point_null").then(|| rp.load(tag));
                    rp.kkt(&mut chk, &e, &mut kkt, a.as_mut());
                    reach.kkt_diag_solves += e.has("a") as usize;
                }
                "starting_basis" => {
                    let mut a = rp.load("initial_ipm");
                    let mut b = Basis::new(control.clone(), model.clone());
                    info = ipx_info(&e.obj("info_in"));
                    starting_basis(&mut a, &mut b, &mut info);
                    check_ipx_info(&mut chk, "info", &e.obj("info"), &info);
                    check_ipx_basis(&mut chk, "basis", &e.obj("basis"), &b);
                    check_ipx_iterate(&mut chk, "iterate", &e.obj("iterate"), &a);
                    reach.forrest_tomlin += (e.obj("basis").str("lu") == "forrest_tomlin") as usize;
                    let si = e.obj("info");
                    reach.dependent += (si.int("dependent_rows") + si.int("dependent_cols")) as i64;
                    it = Some(a);
                    basis = Some(b);
                }
                "basis_ops" => {
                    let b = if e.str("tag") == "starting_basis" {
                        basis.as_ref()
                    } else {
                        basis2.as_ref()
                    };
                    basis_ops(&mut chk, &e, b.expect("basis"));
                }
                "kkt_basis" => {
                    let starting = e.str("tag") == "starting_basis";
                    if !starting {
                        it = Some(rp.load("main_ipm"));
                    }
                    let b = if starting {
                        basis.as_mut()
                    } else {
                        basis2.as_mut()
                    };
                    let mut kkt = KKTSolverBasis::new(control.clone(), b.expect("basis"));
                    info = rp.kkt(&mut chk, &e, &mut kkt, it.as_mut());
                    reach.kkt_basis_solves += e.has("a") as usize;
                    reach.basis_changes += e.int("basis_changes") as i64;
                    let fi = e.obj("info_factorize");
                    reach.dropped += (fi.int("primal_dropped") + fi.int("dual_dropped")) as i64;
                }
                "basis_load" => {
                    let mut b = Basis::new(control.clone(), model.clone());
                    let lerr = b.load(&e.ints("status"));
                    chk.int("errflag", e.int("errflag"), lerr);
                    check_ipx_basis(&mut chk, "basis", &e.obj("basis"), &b);
                    basis2 = Some(b);
                }
                "postsolve" => {
                    postsolve(
                        &mut chk,
                        &e,
                        &model,
                        it.as_mut().expect("iterate"),
                        basis2.as_ref().expect("basis"),
                        &info,
                    );
                    reach.postsolves += 1;
                }
                t => chk.push("t", "known record".into(), t.into(), String::new()),
            }
            Some(chk)
        });
    }
}

/// See ipx_model.rs: `dualize=1` replays in release builds only.
fn skip_in_debug(fx: &Fixture) -> bool {
    cfg!(debug_assertions) && fx.first("params").int("dualize") > 0
}

#[test]
fn ipx_linalg_matches_highs() {
    let mut fails = Failures::default();
    let mut reach = Reach::default();
    for fx in fixtures_in("ipx/linalg") {
        if skip_in_debug(&fx) {
            eprintln!("SKIP {} (dualize=1 in a debug build)", fx.file);
            continue;
        }
        replay(&mut fails, &fx, &mut reach);
    }
    assert!(reach.kkt_diag_solves > 0, "no KKTSolverDiag solve");
    assert!(reach.kkt_basis_solves > 0, "no KKTSolverBasis solve");
    assert!(reach.basis_changes > 0, "no maxvolume basis change");
    assert!(reach.forrest_tomlin > 0, "no ForrestTomlin basis");
    assert!(reach.postsolves > 0, "no postsolve");
    assert!(
        reach.dependent > 0,
        "no dependent rows/cols in a starting basis"
    );
    if !cfg!(debug_assertions) {
        assert!(reach.dualized > 0, "no dualized model");
    }
    eprintln!(
        "ipx_linalg reach: {} kkt_diag solves, {} kkt_basis solves, {} basis changes, \
         {} dropped, {} ForrestTomlin, {} dependent, {} dualized, {} postsolves",
        reach.kkt_diag_solves,
        reach.kkt_basis_solves,
        reach.basis_changes,
        reach.dropped,
        reach.forrest_tomlin,
        reach.dependent,
        reach.dualized,
        reach.postsolves
    );
    fails.finish("ipx_linalg");
}

/// Oracle (b), `ipx/trace/<lp>.<cfg>.jsonl` (`highs-oracle ipx-trace`,
/// highs-oracle-ipx-hooks.patch): the per-iteration IPX trace of a full
/// `Highs::run` with `solver=ipx`. Here only the framing is checked (a
/// regenerated set that lost the hook records fails); the byte diff against
/// the Rust writer is `ipx_run.rs::ipx_trace_matches_highs`.
#[test]
fn ipx_trace_fixtures_are_well_formed() {
    let mut traces = 0;
    for fx in fixtures_in("ipx/trace") {
        assert_eq!(fx.rec(0).t(), "trace_header", "{}", fx.file);
        assert_eq!(fx.rec(fx.len() - 1).t(), "result", "{}", fx.file);
        let ipm: Vec<Rec> = fx.of_type("ipm").collect();
        assert!(!ipm.is_empty(), "{}: no ipm records", fx.file);
        assert_eq!(ipm[0].str("kind"), "start", "{}", fx.file);
        let info = fx.of_type("info").last().expect("info record");
        assert_eq!(info.str("kind"), "solve_end", "{}", fx.file);
        let starts = ipm.iter().filter(|r| r.str("kind") == "start").count();
        let ends = fx
            .of_type("info")
            .filter(|r| r.str("kind") == "solve_end")
            .count();
        assert_eq!((starts, ends), (1, 1), "{}: one LpSolver::Solve", fx.file);
        for (k, r) in fx.recs().skip(1).take(fx.len() - 2).enumerate() {
            assert_eq!(r.int("call"), k as i32, "{}: call sequence", r.ctx);
        }
        traces += 1;
    }
    assert!(traces > 0, "no ipx traces");
}
