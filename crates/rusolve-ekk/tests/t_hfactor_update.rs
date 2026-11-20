//! Golden test for manifest row 24, `src/util/hfactor/update.rs`
//! (HFactor.cpp: updateFT, updatePF, updateMPF, updateAPF, getInvert),
//! entered through `HFactor::update` (row 20), with the post-update
//! FTRAN/BTRAN variants of row 23.
//!
//! Fixtures: `update_sequence_start` / `update_step`* /
//! `update_sequence_end` groups of `factor/<lp>.jsonl` (afiro committed;
//! others target-only), written by `updateSequence` in
//! `tools/highs-oracle/src/factor.cpp`, for update methods FT, PF, MPF, APF
//! from the `slack` and `mid` bases (up to 24 updates, stopping early when
//! `hint` asks for reinversion).
//!
//! Per step k the test checks, in order:
//! 1. `aq`: `aq.clear(); aq.packFlag = true; a.collectAj(aq, variable_in,
//!    1); ftranCall(aq, aq_density)` on the current factor;
//! 2. `ep`: `ep = e_{row_out}; ep.packFlag = true; btranCall(ep,
//!    ep_density)`;
//! 3. `update(&aq, &ep, &row, &hint)` with `row = row_out`, `hint = 0`,
//!    applied to aq/ep rebuilt FROM THE FIXTURE (so a wrong ftran in (1)
//!    does not also corrupt the update), then `hint`, `row_after`,
//!    `aq_after`, `ep_after`;
//! 4. `basic_index[row_out] = variable_in`, then the 4 solves of the step.
//!
//! The factor carries state from step to step, so only the FIRST failing
//! step of a sequence is reported (with every mismatching field of that
//! step); later steps of that sequence are skipped. At the end of a
//! sequence `update_sequence_end.basic_index` and `getInvert()` are
//! compared.

mod common;
mod port;

use common::*;
use port::*;
use rusolve_ekk::prelude::*;

struct Seq<'a> {
    start: Rec<'a>,
    steps: Vec<Rec<'a>>,
    end: Option<Rec<'a>>,
}

fn sequences(fx: &Fixture) -> Vec<Seq<'_>> {
    let mut out: Vec<Seq> = Vec::new();
    for rec in fx.recs() {
        let t = rec.t();
        if t == "update_sequence_start" {
            out.push(Seq {
                start: rec,
                steps: Vec::new(),
                end: None,
            });
            continue;
        }
        if t != "update_step" && t != "update_sequence_end" {
            continue;
        }
        let ctx = rec.ctx.clone();
        let seq = out
            .last_mut()
            .unwrap_or_else(|| panic!("{ctx}: {t} before any update_sequence_start"));
        if t == "update_step" {
            seq.steps.push(rec);
        } else {
            seq.end = Some(rec);
        }
    }
    assert!(!out.is_empty(), "{}: no update sequences", fx.file);
    out
}

fn unit_ep(m: HighsInt, p: HighsInt) -> HVector {
    let mut ep = HVector::default();
    ep.setup(m);
    ep.clear();
    ep.count = 1;
    ep.index[0] = p;
    ep.array[p as usize] = 1.0;
    ep.pack_flag = true;
    ep
}

/// Run one sequence; returns the first failing checker, if any.
fn run_sequence(a: &HighsSparseMatrix, name: &str, seq: &Seq) -> Option<Checker> {
    let s = &seq.start;
    let basis = s.str("basis");
    let method_name = s.str("method");
    let method = s.int("update_method");
    assert_eq!(method, update_method_from_name(method_name), "{}", s.ctx);
    let tag = format!("{name} basis={basis} method={method_name}");
    let m = a.num_row_;

    let mut f = HFactor::default();
    let mut basic_index = s.ints("basic_index");
    setup_factor(&mut f, a, method);
    at("initial HFactor::build".into());
    let rd = f.build(&mut basic_index, None);
    let mut chk = Checker::new(format!(
        "{} {tag} initial build (HFactor::build; see t_hfactor_build)",
        s.ctx
    ));
    chk.int("rank_deficiency", s.int("rank_deficiency"), rd);
    chk.ints(
        "basic_index_after",
        &s.ints("basic_index_after"),
        &basic_index,
    );
    if !chk.ok() {
        return Some(chk);
    }
    if rd != 0 {
        return None;
    }

    for st in &seq.steps {
        let k = st.int("k");
        let q = st.int("variable_in");
        let p = st.int("row_out");
        let mut chk = Checker::new(format!(
            "{} {tag} step k={k} (variable_in={q}, row_out={p}; steps < k passed)",
            st.ctx
        ));
        chk.int(
            "variable_out (= test's basic_index[row_out])",
            st.int("variable_out"),
            basic_index[p as usize],
        );

        at(format!("step k={k} (1) collectAj + HFactor::ftranCall"));
        // (1) aq as HEkkDual::updateFtran forms it.
        let mut aq = HVector::default();
        aq.setup(m);
        aq.clear();
        aq.pack_flag = true;
        a.collect_aj(&mut aq, q, 1.0);
        f.ftran_call(&mut aq, st.f64("aq_density"), None);
        check_hv(
            &mut chk,
            "aq [collectAj + HFactor::ftranCall]",
            &st.hv("aq"),
            &aq,
        );

        at(format!("step k={k} (2) HFactor::btranCall"));
        // (2) ep = B^-T e_p.
        let mut ep = unit_ep(m, p);
        f.btran_call(&mut ep, st.f64("ep_density"), None);
        check_hv(&mut chk, "ep [HFactor::btranCall]", &st.hv("ep"), &ep);

        at(format!(
            "step k={k} (3) HFactor::update -> update{method_name}"
        ));
        // (3) update on the oracle's aq / ep.
        let mut aq_fx = hvector_from_rec(&st.hv("aq"));
        let mut ep_fx = hvector_from_rec(&st.hv("ep"));
        let mut row = p;
        let mut hint: HighsInt = 0;
        f.update(
            &mut basic_index,
            &mut aq_fx,
            &mut ep_fx,
            &mut row,
            &mut hint,
        );
        let upd = format!("[HFactor::update -> update{method_name}]");
        chk.int(&format!("hint {upd}"), st.int("hint"), hint);
        chk.int(&format!("row_after {upd}"), st.int("row_after"), row);
        check_hv(
            &mut chk,
            &format!("aq_after {upd}"),
            &st.hv("aq_after"),
            &aq_fx,
        );
        check_hv(
            &mut chk,
            &format!("ep_after {upd}"),
            &st.hv("ep_after"),
            &ep_fx,
        );

        at(format!("step k={k} (4) post-update solves"));
        // (4) basis change, then the post-update solves.
        basic_index[p as usize] = q;
        check_solves(&mut chk, "solves", &st.solves("solves"), &mut f);
        if !chk.ok() {
            return Some(chk);
        }
        if hint != 0 {
            break;
        }
    }

    let e = seq
        .end
        .as_ref()
        .unwrap_or_else(|| panic!("{}: sequence has no update_sequence_end", s.ctx));
    let mut chk = Checker::new(format!(
        "{} {tag} after {} update(s) (HFactor::getInvert)",
        e.ctx,
        seq.steps.len()
    ));
    chk.ints("basic_index", &e.ints("basic_index"), &basic_index);
    check_invert(&mut chk, "invert", &e.invert("invert"), &f.get_invert());
    (!chk.ok()).then_some(chk)
}

fn run_method(method_name: &str) {
    let mut fails = Failures::default();
    for fx in fixtures_in("factor") {
        let (a, name) = factor_matrix(&fx);
        for seq in sequences(&fx) {
            if seq.start.str("method") != method_name {
                continue;
            }
            let ctx = format!(
                "{} {name} basis={} method={method_name}",
                seq.start.ctx,
                seq.start.str("basis")
            );
            fails.run(&ctx, || run_sequence(&a, &name, &seq));
        }
    }
    fails.finish(&format!("t_hfactor_update method={method_name}"));
}

#[test]
fn update_ft_sequences_match_highs() {
    run_method("FT");
}

// HiGHS guards the PF / MPF / APF branches of ftranL/btranL/ftranU/btranU
// with `assert(update_method != kUpdateMethod...)`; the oracle was built
// with -DNDEBUG (CMAKE_BUILD_TYPE=Release), so those asserts were off when
// the fixtures were written. The port translates them as `debug_assert!`,
// so these three targets only replay HiGHS in a build without debug
// assertions: `cargo test -p rusolve-ekk --release --test t_hfactor_update`.

#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "HiGHS assert(update_method != kUpdateMethodPf); oracle built with NDEBUG: run with --release"
)]
fn update_pf_sequences_match_highs() {
    run_method("PF");
}

#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "HiGHS assert(update_method != kUpdateMethodMpf); oracle built with NDEBUG: run with --release"
)]
fn update_mpf_sequences_match_highs() {
    run_method("MPF");
}

#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "HiGHS assert(update_method != kUpdateMethodApf); oracle built with NDEBUG: run with --release"
)]
fn update_apf_sequences_match_highs() {
    run_method("APF");
}
