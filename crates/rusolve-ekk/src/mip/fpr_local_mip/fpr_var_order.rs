//! Translated from mip-heuristics (MIT License, Copyright (c) 2026 Simon Spoorendonk), commit 36f7d9d. See THIRD_PARTY_NOTICES.md.
//!
//! variable ordering strategies for FPR (fpr_var_order.cpp/h).

use super::VarStrategy;
use crate::HighsInt;
use crate::mip::fpr_local_mip::clique_cover::{build_clique_cover, cliques2_order};
use crate::mip::fpr_local_mip::rng::{S5Rng, s5_shuffle, s5_std_sort, s5_uniform_real};
use crate::mip::mip_solver::HighsMipSolver;
use crate::util::highs_utils::{cpp_max, cpp_min};

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::bucket_by_type (lines 41-56)
fn bucket_by_type(ms: &HighsMipSolver) -> (Vec<HighsInt>, Vec<HighsInt>, Vec<HighsInt>) {
    let model = crate::mip_model!(ms);
    let ncol = model.num_col_;
    let mut bin = Vec::new();
    let mut gen_int = Vec::new();
    let mut cont = Vec::new();

    for j in 0..ncol {
        if !crate::mip::fpr_local_mip::rng::is_integer(&model.integrality_, j) {
            cont.push(j);
        } else if ms.get_domain_2().is_binary(ms, j) {
            bin.push(j);
        } else {
            gen_int.push(j);
        }
    }
    (bin, gen_int, cont)
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::concat_buckets (lines 58-65)
fn concat_buckets(bin: &[HighsInt], gen_int: &[HighsInt], cont: &[HighsInt]) -> Vec<HighsInt> {
    let mut order = Vec::with_capacity(bin.len() + gen_int.len() + cont.len());
    order.extend_from_slice(bin);
    order.extend_from_slice(gen_int);
    order.extend_from_slice(cont);
    order
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::rank_lr (lines 68-73)
fn rank_lr(ms: &HighsMipSolver) -> Vec<HighsInt> {
    let ncol = crate::mip_model!(ms).num_col_;
    (0..ncol).collect()
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::rank_type (lines 76-79)
fn rank_type(ms: &HighsMipSolver) -> Vec<HighsInt> {
    let (bin, gen_int, cont) = bucket_by_type(ms);
    concat_buckets(&bin, &gen_int, &cont)
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::rank_random (lines 82-88)
fn rank_random(ms: &HighsMipSolver, rng: &mut S5Rng) -> Vec<HighsInt> {
    let (mut bin, mut gen_int, mut cont) = bucket_by_type(ms);
    s5_shuffle(&mut bin, rng);
    s5_shuffle(&mut gen_int, rng);
    s5_shuffle(&mut cont, rng);
    concat_buckets(&bin, &gen_int, &cont)
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::rank_locks (lines 91-106)
fn rank_locks(ms: &HighsMipSolver) -> Vec<HighsInt> {
    let up = &ms.mipdata_.uplocks;
    let dn = &ms.mipdata_.downlocks;

    let (mut bin, mut gen_int, mut cont) = bucket_by_type(ms);

    // Stable sort preserves formulation order as tiebreak
    // Comparator: max(up[a], dn[a]) > max(up[b], dn[b])
    let key = |x: HighsInt| -> HighsInt { cpp_max(up[x as usize], dn[x as usize]) };

    bin.sort_by_key(|&x| std::cmp::Reverse(key(x)));
    gen_int.sort_by_key(|&x| std::cmp::Reverse(key(x)));
    cont.sort_by_key(|&x| std::cmp::Reverse(key(x)));

    concat_buckets(&bin, &gen_int, &cont)
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::rank_typecl (lines 122-136)
fn rank_typecl(ms: &HighsMipSolver) -> Vec<HighsInt> {
    let (mut bin, gen_int, cont) = bucket_by_type(ms);
    if bin.is_empty() {
        return concat_buckets(&bin, &gen_int, &cont);
    }

    let clq = &ms.mipdata_.cliquetable;
    let ncol = crate::mip_model!(ms).num_col_;
    let cover = build_clique_cover(&clq.cliques, &clq.cliqueentries, &bin, ncol);

    bin = cover.members;
    bin.extend_from_slice(&cover.uncovered);
    concat_buckets(&bin, &gen_int, &cont)
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::rank_cliques (lines 156-194)
fn rank_cliques(ms: &HighsMipSolver, rng: &mut S5Rng, lp_ref: Option<&[f64]>) -> Vec<HighsInt> {
    if lp_ref.is_none() {
        return rank_typecl(ms);
    }
    let lp_ref = lp_ref.unwrap();

    let (mut bin, gen_int, cont) = bucket_by_type(ms);
    if bin.is_empty() {
        return concat_buckets(&bin, &gen_int, &cont);
    }

    let clq = &ms.mipdata_.cliquetable;
    let ncol = crate::mip_model!(ms).num_col_;
    let cover = build_clique_cover(&clq.cliques, &clq.cliqueentries, &bin, ncol);

    bin.clear();
    let mut keyed: Vec<(HighsInt, f64)> = Vec::new();

    for g in 0..cover.num_groups() {
        keyed.clear();
        let start = cover.group_start[g as usize] as usize;
        let end = cover.group_start[(g + 1) as usize] as usize;

        for i in start..end {
            let col = cover.members[i];
            let mut w = if cover.member_pos[i] != 0 {
                lp_ref[col as usize]
            } else {
                1.0 - lp_ref[col as usize]
            };
            w = cpp_max(w, 1e-10); // avoid zero weights
            let u = s5_uniform_real(rng, 1e-15, 1.0);
            let key = u.ln() / w;
            keyed.push((col, key));
        }

        // std::ranges::sort with descending key comparator
        // s5_std_sort expects "less" comparator: a < b
        // We want descending: a > b, so we pass b < a
        s5_std_sort(&mut keyed, |a, b| b.1 < a.1);

        for (col, _key) in &keyed {
            bin.push(*col);
        }
    }
    bin.extend_from_slice(&cover.uncovered);
    concat_buckets(&bin, &gen_int, &cont)
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::rank_cliques2 (lines 210-226)
fn rank_cliques2(ms: &HighsMipSolver, lp_ref: Option<&[f64]>) -> Vec<HighsInt> {
    if lp_ref.is_none() {
        return rank_typecl(ms);
    }
    let lp_ref = lp_ref.unwrap();

    let (mut bin, gen_int, cont) = bucket_by_type(ms);
    if bin.is_empty() {
        return concat_buckets(&bin, &gen_int, &cont);
    }

    let clq = &ms.mipdata_.cliquetable;
    let dom = ms.get_domain_2();
    let ncol = crate::mip_model!(ms).num_col_;

    bin = cliques2_order(
        &clq.cliques,
        &clq.cliqueentries,
        &bin,
        ncol,
        lp_ref,
        &dom.col_lower_,
        &dom.col_upper_,
    );
    concat_buckets(&bin, &gen_int, &cont)
}

/// mip-heuristics (36f7d9d): fpr_var_order.cpp::compute_var_order (lines 230-253)
pub(crate) fn compute_var_order(
    ms: &HighsMipSolver,
    strategy: VarStrategy,
    rng: &mut S5Rng,
    lp_ref: Option<&[f64]>,
) -> Vec<HighsInt> {
    match strategy {
        VarStrategy::kLR => rank_lr(ms),
        VarStrategy::kType => rank_type(ms),
        VarStrategy::kRandom => rank_random(ms, rng),
        VarStrategy::kLocks => rank_locks(ms),
        VarStrategy::kTypecl => rank_typecl(ms),
        VarStrategy::kCliques => rank_cliques(ms, rng, lp_ref),
        VarStrategy::kCliques2 => rank_cliques2(ms, lp_ref),
        VarStrategy::kDomainSize => {
            // Dynamic strategy: initial order is type-based; actual selection
            // happens at each DFS node in fpr_core via find_smallest_domain.
            rank_type(ms)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp_data::hconst::{HighsVarType, MatrixFormat, kHighsInf};
    use crate::lp_data::highs_callback::HighsCallback;
    use crate::lp_data::highs_lp::HighsLp;
    use crate::lp_data::highs_options::HighsOptions;
    use crate::lp_data::hstruct::HighsSolution;
    use crate::mip::Slot;
    use crate::mip::domain::HighsDomain;
    use crate::mip::fpr_local_mip::rng::S5Rng;

    fn set_root_domain(ms: &mut HighsMipSolver) {
        let mut dom = HighsDomain::default();
        dom.col_lower_ = vec![0.0, 0.0, 0.0, 0.0];
        dom.col_upper_ = vec![1.0, 1.0, 10.0, 10.0];
        ms.mipdata_.domains = vec![Slot::new(dom)];
    }

    fn make_small_mip() -> HighsMipSolver {
        let mut lp = HighsLp::default();
        lp.num_col_ = 4;
        lp.num_row_ = 1;
        lp.col_cost_ = vec![1.0, 1.0, 1.0, 1.0];
        lp.col_lower_ = vec![0.0, 0.0, 0.0, 0.0];
        lp.col_upper_ = vec![1.0, 1.0, 10.0, 10.0];
        lp.integrality_ = vec![
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kInteger,
            HighsVarType::kContinuous,
        ];
        lp.row_lower_ = vec![-kHighsInf];
        lp.row_upper_ = vec![10.0];
        lp.a_matrix_.format_ = MatrixFormat::kColwise;
        lp.a_matrix_.num_col_ = 4;
        lp.a_matrix_.num_row_ = 1;
        lp.a_matrix_.start_ = vec![0, 1, 2, 3, 4];
        lp.a_matrix_.index_ = vec![0, 0, 0, 0];
        lp.a_matrix_.value_ = vec![1.0, 1.0, 1.0, 1.0];

        let options = HighsOptions::default();
        let mut cb = HighsCallback::default();
        let solution = HighsSolution::default();
        let mut ms = HighsMipSolver::new(&mut cb, &options, &lp, &solution, false, 0);
        set_root_domain(&mut ms);
        ms
    }

    #[test]
    fn s5_varsel_lr_identity() {
        let ms = make_small_mip();
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kLR, &mut rng, None);
        assert_eq!(order, vec![0, 1, 2, 3]);
    }

    #[test]
    fn s5_varsel_type_buckets() {
        let ms = make_small_mip();
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kType, &mut rng, None);
        // Cols 0,1 are binary (lb=0, ub=1, integer). Col 2 is general integer. Col 3 is continuous.
        // Order should be: bin (0,1), gen_int (2), cont (3)
        assert_eq!(order, vec![0, 1, 2, 3]);
    }

    #[test]
    fn s5_varsel_random_permutation() {
        let ms = make_small_mip();
        let mut rng = S5Rng::new(123);
        let order = compute_var_order(&ms, VarStrategy::kRandom, &mut rng, None);
        assert_eq!(order.len(), 4);
        let mut sorted = order.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2, 3]);
    }

    #[test]
    fn s5_varsel_locks_sort() {
        let mut ms = make_small_mip();
        ms.mipdata_.uplocks = vec![1, 5, 2, 0];
        ms.mipdata_.downlocks = vec![2, 1, 3, 0];
        // max(up, dn): col0=2, col1=5, col2=3, col3=0
        // Binary bucket: 0, 1. Sorted by max desc: 1 (5), 0 (2).
        // Gen int: 2. Cont: 3.
        // Result: 1, 0, 2, 3
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kLocks, &mut rng, None);
        assert_eq!(order, vec![1, 0, 2, 3]);
    }

    #[test]
    fn s5_varsel_typecl_no_cliques() {
        let ms = make_small_mip();
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kTypecl, &mut rng, None);
        // No cliques in this simple model, so it should fall back to type order
        assert_eq!(order, vec![0, 1, 2, 3]);
    }

    #[test]
    fn s5_varsel_cliques_none_fallback() {
        let ms = make_small_mip();
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kCliques, &mut rng, None);
        // lp_ref is None, so it falls back to typecl (which falls back to type)
        assert_eq!(order, vec![0, 1, 2, 3]);
    }

    #[test]
    fn s5_varsel_cliques2_none_fallback() {
        let ms = make_small_mip();
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kCliques2, &mut rng, None);
        // lp_ref is None, so it falls back to typecl (which falls back to type)
        assert_eq!(order, vec![0, 1, 2, 3]);
    }

    fn add_clique(ms: &mut HighsMipSolver, cols: &[HighsInt]) {
        use crate::HighsUInt;
        use crate::lp_data::hconst::kHighsIInf;
        use crate::mip::clique_table::{HighsCliqueTable_Clique, HighsCliqueTable_CliqueVar};
        let clq = &mut ms.mipdata_.cliquetable;
        let mut cl = HighsCliqueTable_Clique::default();
        cl.start = clq.cliqueentries.len() as HighsInt;
        cl.end = cl.start + cols.len() as HighsInt;
        cl.origin = kHighsIInf;
        cl.num_zero_fixed = 0;
        cl.equality = false;
        clq.cliques.push(cl);
        for &c in cols {
            clq.cliqueentries
                .push(HighsCliqueTable_CliqueVar::from_parts(c as HighsUInt, 1));
        }
    }

    #[test]
    fn s5_varsel_typecl_uses_clique_entry_order() {
        let mut ms = make_small_mip();
        add_clique(&mut ms, &[1, 0]);
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kTypecl, &mut rng, None);
        assert_eq!(order, vec![1, 0, 2, 3]);
    }

    #[test]
    fn s5_varsel_cliques_efraimidis_spirakis_keys() {
        let mut ms = make_small_mip();
        add_clique(&mut ms, &[1, 0]);
        let lp = vec![0.3, 0.7, 4.0, 2.0];
        let mut rng = S5Rng::new(7);
        let order = compute_var_order(&ms, VarStrategy::kCliques, &mut rng, Some(&lp));

        // Same draws in the same order: group member 1 first, then member 0.
        let mut r2 = S5Rng::new(7);
        let k1 = s5_uniform_real(&mut r2, 1e-15, 1.0).ln() / 0.7;
        let k0 = s5_uniform_real(&mut r2, 1e-15, 1.0).ln() / 0.3;
        let expect_bin = if k0 > k1 { vec![0, 1] } else { vec![1, 0] };
        assert_eq!(&order[..2], expect_bin.as_slice());
        assert_eq!(&order[2..], &[2, 3]);
        // Exactly two draws consumed.
        assert_eq!(rng.next_u64(), r2.next_u64());
    }

    #[test]
    fn s5_varsel_cliques2_most_positive_literal_first() {
        let mut ms = make_small_mip();
        add_clique(&mut ms, &[0, 1]);
        let lp = vec![0.2, 0.8, 4.0, 2.0];
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kCliques2, &mut rng, Some(&lp));
        assert_eq!(order, vec![1, 0, 2, 3]);
    }

    #[test]
    fn s5_varsel_domain_size_type() {
        let ms = make_small_mip();
        let mut rng = S5Rng::new(42);
        let order = compute_var_order(&ms, VarStrategy::kDomainSize, &mut rng, None);
        // kDomainSize initial order is type-based
        assert_eq!(order, vec![0, 1, 2, 3]);
    }
}
