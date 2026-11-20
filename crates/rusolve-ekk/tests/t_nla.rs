//! Golden test for manifest rows 27-28, `src/simplex/nla/{product_form,
//! nla}.rs` (HSimplexNla.cpp: setup, invert, ftran, btran, update,
//! ftranInScaledSpace, btranInScaledSpace; HSimplexNlaProductForm.cpp on the
//! non-frozen path).
//!
//! Fixtures: the HFactor records of `factor/<lp>.jsonl` (afiro committed;
//! others target-only). HEkk factorises the *scaled* LP, and once the LP is
//! scaled in place (`is_scaled_ = true`) HSimplexNla has no scale to apply
//! (`scale_ = NULL`), so `HSimplexNla::setup(lp, basic_index, options, a,
//! kDefaultPivotThreshold)` + `invert()` + `ftran`/`btran` + `update` must
//! reproduce the oracle's bare-HFactor results exactly: setup passes
//! num_row as num_basic, `kDefaultPivotThreshold`, `options.
//! factor_pivot_tolerance` (= kDefaultPivotTolerance), `highs_debug_level`
//! (= kHighsDebugLevelMin), use_original_HFactor_logic = true and
//! kUpdateMethodFt, i.e. the oracle's `setupFactor(..., kUpdateMethodFt)`.
//!
//! Compared: every `build` record (via `invert`, plus
//! `build_synthetic_tick_`, the MUST-PORT tick copy), its HVector solves,
//! and the FT `update_step` sequences (HFactor-level details as in
//! `t_hfactor_update`).
//!
//! Not covered (fixture gap G-NLA-SCALE): the scaled-space transforms with
//! `scale_ != NULL` (`applyBasisMatrixRowScale`, `transformForUpdate`,
//! `basicColScaleFactor`, `pivotInScaledSpace`, `rowEp2NormInScaledSpace`)
//! and the frozen-basis product form; they are reached in the simplex only
//! through HEkk and are gated by `t_ekk_init` / `t_dual_trace`.

mod common;
mod port;

use common::*;
use port::*;
use rusolve_ekk::prelude::*;

struct Nla {
    nla: HSimplexNla,
    lp: HighsLp,
    options: HighsOptions,
    basic_index: Vec<HighsInt>,
}

impl Nla {
    fn new(a: &HighsSparseMatrix, basic_index: &[HighsInt]) -> Nla {
        let lp = HighsLp {
            num_col_: a.num_col_,
            num_row_: a.num_row_,
            a_matrix_: a.clone(),
            is_scaled_: true,
            ..HighsLp::default()
        };
        let options = HighsOptions::default();
        let mut nla = HSimplexNla::default();
        nla.setup(&lp, basic_index, &options, a, kDefaultPivotThreshold);
        Nla {
            nla,
            lp,
            options,
            basic_index: basic_index.to_vec(),
        }
    }
}

impl Solver for Nla {
    fn name(&self) -> &'static str {
        "HSimplexNla::ftran/btran (scale_ = NULL)"
    }
    fn hv(&mut self, ftran: bool, v: &mut HVector, ed: f64) {
        let refs = NlaRefs::new(&self.lp, &self.basic_index, &self.options);
        assert!(
            refs.scale_.is_none(),
            "is_scaled_ LP must give scale_ = NULL"
        );
        if ftran {
            self.nla.ftran(&refs, v, ed, None);
        } else {
            self.nla.btran(&refs, v, ed, None);
        }
    }
    fn vec(&mut self, _ftran: bool, _v: &mut Vec<f64>) {
        unreachable!("HSimplexNla has no std::vector overload; callers skip *_vec solves")
    }
}

fn hv_solves(rec: &Rec) -> Vec<SolveRec> {
    rec.solves("solves")
        .into_iter()
        .filter(|s| !s.op.ends_with("_vec"))
        .collect()
}

#[test]
fn nla_invert_and_solves_match_highs() {
    let mut fails = Failures::default();
    let options_check = HighsOptions::default();
    assert_eq!(
        options_check.factor_pivot_tolerance.to_bits(),
        kDefaultPivotTolerance.to_bits(),
        "HighsOptions::default().factor_pivot_tolerance must be kDefaultPivotTolerance"
    );
    for fx in fixtures_in("factor") {
        let (a, name) = factor_matrix(&fx);
        for rec in recs_of(&fx, "build") {
            let basis = rec.str("basis");
            let ctx = format!(
                "{} {name} basis={basis} (HSimplexNla::setup + invert)",
                rec.ctx
            );
            fails.run(&ctx.clone(), || {
                let mut s = Nla::new(&a, &rec.ints("basic_index"));
                let mut bi = s.basic_index.clone();
                at("HSimplexNla::invert".into());
                let rd = s.nla.invert(&mut bi);
                s.basic_index = bi.clone();
                let mut chk = Checker::new(ctx.clone());
                let e = rec.build_result("result");
                check_build_result(&mut chk, &e, &s.nla.factor_, rd, &bi);
                chk.f64(
                    "HSimplexNla::build_synthetic_tick_ (= factor_.build_synthetic_tick)",
                    e.build_synthetic_tick,
                    s.nla.build_synthetic_tick_,
                );
                at("HSimplexNla::ftran/btran".into());
                check_solves(&mut chk, "solves", &hv_solves(&rec), &mut s);
                (!chk.ok()).then_some(chk)
            });
        }
    }
    fails.finish("t_nla invert + solves");
}

#[test]
fn nla_ft_update_sequences_match_highs() {
    let mut fails = Failures::default();
    for fx in fixtures_in("factor") {
        let (a, name) = factor_matrix(&fx);
        let recs: Vec<Rec> = fx.recs().collect();
        for (i, start) in recs.iter().enumerate() {
            if start.t() != "update_sequence_start" || start.str("method") != "FT" {
                continue;
            }
            let basis = start.str("basis");
            let ctx = format!("{} {name} basis={basis} method=FT (HSimplexNla)", start.ctx);
            let steps: Vec<&Rec> = recs[i + 1..]
                .iter()
                .take_while(|r| r.t() == "update_step")
                .collect();
            let end = &recs[i + 1 + steps.len()];
            assert_eq!(end.t(), "update_sequence_end", "{}", end.ctx);
            fails.run(&ctx.clone(), || {
                let mut s = Nla::new(&a, &start.ints("basic_index"));
                let mut bi = s.basic_index.clone();
                at("HSimplexNla::invert".into());
                let rd = s.nla.invert(&mut bi);
                s.basic_index = bi;
                let mut chk = Checker::new(format!("{ctx} initial invert"));
                chk.int("rank_deficiency", start.int("rank_deficiency"), rd);
                chk.ints(
                    "basic_index_after",
                    &start.ints("basic_index_after"),
                    &s.basic_index,
                );
                if !chk.ok() || rd != 0 {
                    return (!chk.ok()).then_some(chk);
                }
                let m = a.num_row_;
                for st in &steps {
                    let k = st.int("k");
                    let q = st.int("variable_in");
                    let p = st.int("row_out");
                    let mut chk = Checker::new(format!(
                        "{} {name} basis={basis} method=FT step k={k} (variable_in={q}, row_out={p}; steps < k passed)",
                        st.ctx
                    ));
                    at(format!("step k={k} aq: collectAj + HSimplexNla::ftran"));
                    let mut aq = HVector::default();
                    aq.setup(m);
                    aq.clear();
                    aq.pack_flag = true;
                    a.collect_aj(&mut aq, q, 1.0);
                    s.hv(true, &mut aq, st.f64("aq_density"));
                    check_hv(&mut chk, "aq [HSimplexNla::ftran]", &st.hv("aq"), &aq);

                    at(format!("step k={k} ep: HSimplexNla::btran"));
                    let mut ep = HVector::default();
                    ep.setup(m);
                    ep.clear();
                    ep.count = 1;
                    ep.index[0] = p;
                    ep.array[p as usize] = 1.0;
                    ep.pack_flag = true;
                    s.hv(false, &mut ep, st.f64("ep_density"));
                    check_hv(&mut chk, "ep [HSimplexNla::btran]", &st.hv("ep"), &ep);

                    at(format!("step k={k} HSimplexNla::update"));
                    let mut aq_fx = hvector_from_rec(&st.hv("aq"));
                    let mut ep_fx = hvector_from_rec(&st.hv("ep"));
                    let mut row = p;
                    let mut hint: HighsInt = 0;
                    let mut bi = s.basic_index.clone();
                    s.nla
                        .update(&mut bi, &mut aq_fx, &mut ep_fx, &mut row, &mut hint);
                    s.basic_index = bi;
                    chk.int("hint [HSimplexNla::update]", st.int("hint"), hint);
                    chk.int("row_after [HSimplexNla::update]", st.int("row_after"), row);
                    check_hv(
                        &mut chk,
                        "aq_after [HSimplexNla::update]",
                        &st.hv("aq_after"),
                        &aq_fx,
                    );
                    check_hv(
                        &mut chk,
                        "ep_after [HSimplexNla::update]",
                        &st.hv("ep_after"),
                        &ep_fx,
                    );
                    s.basic_index[p as usize] = q;
                    at(format!("step k={k} post-update HSimplexNla::ftran/btran"));
                    check_solves(&mut chk, "solves", &hv_solves(st), &mut s);
                    if !chk.ok() {
                        return Some(chk);
                    }
                    if hint != 0 {
                        break;
                    }
                }
                let mut chk = Checker::new(format!("{} {name} basis={basis} method=FT end", end.ctx));
                chk.ints("basic_index", &end.ints("basic_index"), &s.basic_index);
                check_invert(&mut chk, "invert", &end.invert("invert"), &s.nla.factor_.get_invert());
                (!chk.ok()).then_some(chk)
            });
        }
    }
    fails.finish("t_nla FT update sequences");
}
