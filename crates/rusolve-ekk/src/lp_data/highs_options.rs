//! HiGHS source: lp_data/HighsOptions.h (subset)
//!
//! Translated from HiGHS (MIT License, Copyright (c) 2026 HiGHS). See THIRD_PARTY_NOTICES.md.
//!
//! Only the options read by the simplex path (HEkk*, HSimplex*, HFactor,
//! HighsSparseMatrix), presolve, the IPX wrapper, the MIP solver (`mip/`) and PDLP (`pdlp/`). Every default is the HiGHS default from the matching
//! `OptionRecord*` in `HighsOptions.h` (the third numeric argument of an
//! `OptionRecordInt/Double(name, ..., &var, lower, default, upper)` and the
//! last argument of `OptionRecordBool/String`). The option-record machinery
//! (validation, setters by name, file I/O) is not ported.

use crate::HighsInt;
use crate::io::highs_io::HighsLogOptions;
use crate::lp_data::hconst::*;
use crate::simplex::simplex_const::*;

/// HiGHS: HighsOptions.h::HighsOptionsStruct / HighsOptions (simplex subset)
#[derive(Clone, Debug)]
pub struct HighsOptions {
    // Run-time / general
    pub time_limit: f64,
    pub objective_bound: f64,
    pub objective_target: f64,
    pub threads: HighsInt,
    pub random_seed: HighsInt,
    pub parallel: String,
    pub highs_debug_level: HighsInt,
    pub highs_analysis_level: HighsInt,
    pub allow_unbounded_or_infeasible: bool,

    // Presolve / solver choice (HighsOptions.h, presolve)
    /// "off", "choose" or "on" (HiGHS default "choose").
    pub presolve: String,
    pub solver: String,
    pub run_crossover: String,
    pub use_implied_bounds_from_presolve: bool,
    pub lp_presolve_requires_basis_postsolve: bool,
    pub presolve_reduction_limit: HighsInt,
    pub restart_presolve_reduction_limit: HighsInt,
    pub presolve_substitution_maxfillin: HighsInt,
    pub presolve_rule_off: HighsInt,
    pub presolve_rule_test: HighsInt,
    pub presolve_rule_logging: bool,
    pub presolve_remove_slacks: bool,
    pub presolve_pivot_threshold: f64,
    pub mip_feasibility_tolerance: f64,

    // MIP (HighsOptions.h, MIP; HiGHS 1.15.1 defaults)
    pub mip_detect_symmetry: bool,
    pub mip_allow_restart: bool,
    pub mip_max_nodes: HighsInt,
    pub mip_max_stall_nodes: HighsInt,
    pub mip_max_start_nodes: HighsInt,
    pub mip_debug_solution_file: String,
    pub mip_improving_solution_save: bool,
    pub mip_improving_solution_report_sparse: bool,
    pub mip_improving_solution_file: String,
    pub mip_root_presolve_only: bool,
    pub mip_lifting_for_probing: HighsInt,
    pub mip_max_leaves: HighsInt,
    pub mip_max_improving_sols: HighsInt,
    pub mip_lp_age_limit: HighsInt,
    pub mip_pool_age_limit: HighsInt,
    pub mip_pool_soft_limit: HighsInt,
    pub mip_pscost_minreliable: HighsInt,
    pub mip_min_cliquetable_entries_for_parallelism: HighsInt,
    pub mip_report_level: HighsInt,
    pub mip_heuristic_effort: f64,
    pub mip_heuristic_run_feasibility_jump: bool,
    pub mip_heuristic_run_rins: bool,
    pub mip_heuristic_run_rens: bool,
    pub mip_heuristic_run_root_reduced_cost: bool,
    pub mip_heuristic_run_zi_round: bool,
    pub mip_heuristic_run_shifting: bool,
    pub mip_allow_cut_separation_at_nodes: bool,
    pub mip_rel_gap: f64,
    pub mip_abs_gap: f64,
    pub mip_min_logging_interval: f64,
    pub mip_lp_solver: String,
    pub mip_ipm_solver: String,
    pub mip_search_simulate_concurrency: bool,
    // rusolve backport: HiGHS latest #3213, not a 1.15.1 option
    pub mip_sb_budget_pre_restart: bool,
    // rusolve backport: HiGHS latest #3445, not a 1.15.1 option
    pub mip_dfs_node_limit: bool,
    // rusolve backport: HiGHS latest #3336, not a 1.15.1 option
    pub mip_reduced_stack_fast: bool,
    // rusolve backport: HiGHS latest #3184, not a 1.15.1 option
    pub mip_cdouble_act_contrib: bool,
    // rusolve backport: HiGHS latest #3357, not a 1.15.1 option
    pub mip_keep_branch_reason: bool,
    // rusolve backport: HiGHS latest #3177, not a 1.15.1 option
    pub mip_redundant_vbd: bool,
    // rusolve backport: HiGHS latest #3179, not a 1.15.1 option
    pub mip_path_same_bds: bool,
    // rusolve backport: HiGHS latest #3468, not a 1.15.1 option
    pub mip_strengthen_small_contrib: bool,
    // rusolve backport: HiGHS latest #3474/#3481, not a 1.15.1 option
    pub mip_lurking_bounds_robust: bool,
    // rusolve backport: HiGHS latest #3404, not a 1.15.1 option
    pub mip_dual_subst_limits: bool,
    // rusolve backport: HiGHS latest #3437, not a 1.15.1 option
    pub mip_ranged_dual_implied: bool,
    // rusolve backport: HiGHS latest #2990, not a 1.15.1 option
    pub mip_machine_sched_sepa: bool,
    // rusolve backport: HiGHS latest #3244, not a 1.15.1 option
    pub mip_dual_fix_probing: bool,
    // rusolve backport: HiGHS latest #3264, not a 1.15.1 option
    pub mip_aggregate_implications: bool,
    // rusolve backport: HiGHS latest #3286, not a 1.15.1 option
    pub mip_clique_normalise: bool,
    // rusolve backport: HiGHS 1.16.0 #3318, not a 1.15.1 option
    pub mip_coef_change_impl_bounds: bool,
    // rusolve backport: HiGHS 1.16.0 #3332, not a 1.15.1 option
    pub mip_dual_side_lost: bool,
    // rusolve backport: HiGHS latest #3223, not a 1.15.1 option
    pub mip_cdouble_proofs: bool,
    // rusolve backport: HiGHS latest #3220, not a 1.15.1 option
    pub mip_fix_last_active: bool,
    // rusolve backport: HiGHS latest #3209, not a 1.15.1 option
    pub mip_keep_implications: bool,
    // rusolve backport: HiGHS 1.16.0 #3312, not a 1.15.1 option
    pub mip_presolve_clique_table: bool,
    // rusolve backport: HiGHS 1.16.0 #3218 + #3250 (row/column append and
    // cut-row types in the postsolve stack, moveCutsToPool), not a 1.15.1 option
    pub mip_presolve_row_append: bool,
    // rusolve backport: HiGHS 1.16.0 #3255 Fourier-Motzkin presolve, not a
    // 1.15.1 option; effective only together with `mip_presolve_row_append` and
    // `h16_presolve_rules`
    pub mip_fm_presolve: bool,
    // rusolve backport: HiGHS 1.16.0 #3535 (u64 symmetry vertex
    // hashes, symmetry-detection time limit), not a 1.15.1 option
    pub mip_sym_hash64: bool,
    // rusolve backport: HiGHS 1.16.0 options, not 1.15.1 options; all default false
    pub h16_ac_centring: bool,
    pub h16_ac_fix: bool,
    pub h16_add_no_modify_on_error: bool,
    pub h16_agg_term_zero: bool,
    pub h16_assess_small: bool,
    pub h16_big_act: bool,
    pub h16_bound_chain: bool,
    pub h16_bound_range_inf: bool,
    pub h16_cdouble: bool,
    pub h16_clique_links: bool,
    pub h16_col_stuffing: bool,
    pub h16_cost_perturb: bool,
    pub h16_crossover_x: bool,
    pub h16_cupdlp_time: bool,
    pub h16_cutpool_age: bool,
    pub h16_dep_eq_time_cap: bool,
    pub h16_dominated_cols: bool,
    pub h16_dual_free_list: bool,
    pub h16_dual_postsolve_rules: bool,
    pub h16_dup_row_undo: bool,
    pub h16_dupcol_warn: bool,
    pub h16_enum_sig64: bool,
    pub h16_eq_row_adds: bool,
    pub h16_feasrelax_sense: bool,
    pub h16_fj: bool,
    pub h16_frontier_merge: bool,
    pub h16_hash_common: bool,
    pub h16_hfactor_pivot: bool,
    pub h16_hipdlp_max: bool,
    pub h16_hipo_display: bool,
    pub h16_hipo_empty_rows: bool,
    pub h16_impl_eq: bool,
    pub h16_implint_bounds: bool,
    pub h16_int_bound_round: bool,
    pub h16_int_bound_rounding: bool,
    pub h16_int_round_assert_tol: bool,
    pub h16_ipm_gap_offset: bool,
    pub h16_ipx_free_rows: bool,
    pub h16_ipx_scaled: bool,
    pub h16_ipx_stall: bool,
    pub h16_kkt_effective_costs: bool,
    pub h16_max_duals: bool,
    pub h16_mip_start: bool,
    pub h16_mip_workers: bool,
    pub h16_ninfmin: bool,
    pub h16_nla_debug_gate: bool,
    pub h16_obj_cdouble: bool,
    pub h16_parallel_rows: bool,
    pub h16_postsolve_basis_check: bool,
    pub h16_presolve_rules: bool,
    /// HiGHS 1.16 (#3536): HighsOptions.h::presolve_mode (1.16.0 lines 1046-1050), read
    /// only under `h16_presolve_rules`.
    pub presolve_mode: HighsInt,
    /// HiGHS 1.16 (#3255): HighsOptions.h::presolve_fm_obj_reformulation, read only by the
    /// Fourier-Motzkin presolve (`mip_fm_presolve`).
    pub presolve_fm_obj_reformulation: bool,
    pub h16_row_duals: bool,
    pub h16_fourier_motzkin: bool,
    pub h16_qp_bounded_start: bool,
    pub h16_qp_convexity: bool,
    pub h16_ray_simplex: bool,
    pub h16_rebuild_reason: bool,
    pub h16_repair_lp: bool,
    pub h16_row_presolve: bool,
    pub h16_row_redundant: bool,
    /// HiGHS 1.16 (#3350 and the 1.16.0 coefficient-range report): HighsSolve.cpp::
    /// assessExcessiveObjectiveBoundScaling reports effective cost ranges, split
    /// continuous/non-continuous bound ranges, the 1.16.0 warnings and the
    /// tolerance-ratio warnings. Log only. Off keeps the 1.15.1 report.
    pub h16_scaling_report: bool,
    pub h16_singleton_loops: bool,
    pub h16_subst_int: bool,
    pub h16_sym_refine: bool,
    pub h16_unresolved_unbounded_or_infeasible: bool,
    pub h16_user_scale_offset: bool,
    pub h16_user_sol_obj: bool,
    pub h16_user_solution_assess: bool,
    pub h16_watch_link: bool,
    pub h16_zero_cost_singleton: bool,
    /// rusolve port plumbing (no HiGHS PR): the solvers the MIP
    /// builds internally (the LP relaxation's `lpsolver`, the analytic-centre IPM, the
    /// repair LP and the root IPM fallback) start from `Highs::default()` as in C++, which
    /// in the port means 1.15.1 code. In HiGHS 1.16.0 they run 1.16.0 code. On: each such
    /// solver copies the parent's `h16_*` / `mip_*` / `mip_*` flags ([`HighsOptions::copy_version_flags_from`]).
    /// Off reproduces the 1.15.1-code sub-solvers.
    pub h16_sub_solver_flags: bool,
    // rusolve: SCIP 11 conf_usegenres, not a 1.15.1 option
    pub mip_conflict_resolution: bool,
    /// rusolve, only read when `mip_conflict_resolution` is
    /// on: skip the generalized-resolution step while `mipdata_.num_nodes == 0`, i.e. in the root's
    /// heuristic dives and in every (sub-)MIP before its first node is counted. The 1.15.1
    /// graph-based analysis still runs there. On decomp2 resolution at the root changes the RENS
    /// sub-MIP (116 nodes vs 1); with this gate the trace is the 1.16.0 one on seeds 0-2.
    pub mip_conflict_resolution_tree_only: bool,
    // rusolve: SCIP 11 propagating/symmetry/usedynamicprop, not a 1.15.1 option
    pub mip_symmetry_dynamic: bool,
    // rusolve: SCIP 11 cutselection/hybrid, not a 1.15.1 option
    pub mip_cut_selection_hybrid: bool,
    pub mip_knapsack_cover_sepa: bool,
    pub mip_dff_separator: bool,
    // rusolve: relax the tableau separator's weight-range and density limits
    pub mip_tableau_loosen: bool,
    pub mip_tableau_loosen_auto: bool,
    pub mip_ti_precedence_sepa: bool,
    pub mip_ti_list_sched: bool,
    pub mip_skip_analytic_centre: bool,
    pub mip_netdesign_sepa: bool,
    pub mip_netdesign_heur: bool,
    pub mip_netdesign_cutset: bool,
    // rusolve: restart so presolve tightens
    // big-M coefficients against the objective cutoff; not a HiGHS option
    pub mip_cutoff_bigm: bool,
    pub mip_node_bigm: bool,
    pub mip_node_bigm_rows: bool,
    pub mip_node_bigm_min_ratio: f64,
    pub mip_node_bigm_active: f64,
    /// rusolve: relative range improvement for propagation to accept a continuous bound
    /// (HiGHS hard-codes 0.3; SCIP numerics/boundstreps is 0.05).
    pub mip_cont_prop_improve: f64,
    /// rusolve: 0 = HiGHS count-based pseudocost reliability. > 0: sample-variance reliability
    /// (Hendel 2015).
    pub mip_pscost_variance_tol: f64,
    /// rusolve: RENS skips its sub-MIP when the fixing rate is below this (0 = off).
    pub mip_rens_min_fixing_rate: f64,
    // rusolve: SCIP 11 heuristics/oneopt, not a 1.15.1 option
    pub mip_heuristic_run_oneopt: bool,
    // rusolve: SCIP 11 heuristics/crossover, not a 1.15.1 option
    pub mip_heuristic_run_crossover: bool,
    // rusolve: mip-heuristics (Spoorendonk, MIT) fpr / local_mip presolve heuristics, not a 1.15.1 option
    pub mip_heuristic_run_fpr: bool,
    // rusolve: mip-heuristics (Spoorendonk, MIT) fpr / local_mip presolve heuristics, not a 1.15.1 option
    pub mip_heuristic_run_local_mip: bool,
    /// rusolve, no HiGHS PR: in the 1.16.0
    /// `strengthenInequalities` cover loop (`mip_strengthen_small_contrib`), stop the loop
    /// for the current row when `slackupper < 0`, the invariant HiGHS asserts but compiles
    /// out. Off reproduces HiGHS 1.16.0, which can loop forever. On in the recommended
    /// profile (`RECOMMENDED_MIP_ON`), off in `Default` (1.15.1).
    pub mip_strengthen_cover_guard: bool,
    // rusolve ls-gub-swap: GUB swap / ejection-chain moves in
    // LocalMIP (Benoist et al. 2011), not a HiGHS option; acts only with run_local_mip
    pub mip_ls_gub_swap: bool,
    // rusolve: fix the integers of LS end points and LP the rest, not a 1.15.1 option
    pub mip_ls_int_lp_repair: bool,
    // rusolve: bandit-scheduled LNS (RENS, RINS,
    // mutation, crossover) on the main thread at the SCIP-heuristic hook sites; not a HiGHS option
    pub mip_alns: bool,
    // rusolve: substitute a unit-coefficient integer column out of a short
    // integer equation, turning its non-implied bounds into one new row; not a HiGHS option
    pub mip_presolve_eq_subst: bool,
    // rusolve: merge variable-bound rows over one clique into one row
    pub mip_presolve_clique_varbound: bool,
    pub mip_presolve_xor_gf2: bool,
    pub mip_orbitope_root_fix: bool,
    pub mip_orbitope_full_fix: bool,
    // rusolve: Shift-and-Propagate, a pre-root LP-free primal
    // heuristic (Berthold & Hendel 2015, SCIP heur_shiftandpropagate.c); not a HiGHS option
    pub mip_shift_and_propagate: bool,
    // rusolve: proximity search (Fischetti-Monaci 2014, SCIP heur_proximity.c) and local branching (Fischetti-Lodi 2003, SCIP heur_localbranching.c) improvement sub-MIPs on the main thread; not HiGHS options
    pub mip_proximity: bool,
    // rusolve: pre-root proximity-search chain right after the feasibility jump
    pub mip_improve_early: bool,
    // rusolve: stall-triggered proximity search (incumbent and dual bound both stalled), no other heuristic option needed
    pub mip_primal_stall: bool,
    // rusolve: work budget of that chain, in multiples of the pre-root work
    pub mip_improve_early_work: f64,
    pub mip_guided_dive: bool,
    pub mip_local_branching: bool,
    // rusolve: LP-free conflict-driven DFS with Luby restarts as a pre-root
    // and post-root primal heuristic on pure-integer models; not a HiGHS option
    pub mip_cdcl_dive: bool,
    /// rusolve: stop root separation after 5 rounds with a flat bound once an incumbent exists.
    pub mip_root_bound_stall: bool,
    // rusolve: first-UIP clauses with two watched literals inside the dive
    pub mip_cdcl_uip: bool,
    // rusolve: CDCL root feasibility search when the objective is zero
    pub mip_zero_obj_cdcl: bool,
    // rusolve: solve connected components of a MIP as separate MIPs
    pub mip_decompose: bool,
    pub mip_lattice: bool,
    pub mip_activation_dff: bool,
    pub mip_fj_native: bool,
    // rusolve: no incumbent after the pre-root FJ -> native FJ, then FPR and Local-MIP
    pub mip_first_inc: bool,
    // rusolve: the same fallback after the root separation and root heuristics
    pub mip_first_inc_root: bool,
    pub mip_lattice_kernel: bool,
    /// Stop MIP presolve probing after 10000 probed binaries.
    pub mip_presolve_probing_cap: bool,
    pub mip_presolve_domcol_full: bool,
    /// rusolve: presolve probing stops after this many successive probes without a deletion or
    /// a new clique (SCIP presolving/probing/maxuseless); 0 = off.
    pub mip_presolve_probing_useless: HighsInt,
    /// rusolve ablation lever: bit mask of conflict-analysis sources
    /// to switch off (1 tree propagation, 2 heuristics, 4 infeasible-LP dual ray, 8 bound-exceeding
    /// LP dual proof, 16 degenerate-dual reconvergence, 32 conflict cuts from dual proofs); 0 = HiGHS behaviour.
    pub mip_conflict_off_mask: HighsInt,
    pub mip_redcost_proof_uncapped: bool,
    pub mip_root_stall_scip: bool,
    pub mip_heur_lp_quota: bool,
    // rusolve: k diversified copies of the MIP search race on k threads and share incumbents; 0/1 = off
    pub mip_racing_threads: HighsInt,
    // rusolve: race the root LP between dual simplex and IPM + crossover on two threads; not a HiGHS option
    pub mip_root_lp_concurrent: bool,
    pub solve_relaxation: bool,
    /// HiGHS: HighsOptions.h `user_objective_scale` (record lines 874-877: default 0,
    /// bounds [-kHighsIInf, kHighsIInf])
    pub user_objective_scale: HighsInt,
    /// HiGHS: HighsOptions.h `user_bound_scale` (record lines 879-882: default 0,
    /// bounds [-kHighsIInf, kHighsIInf])
    pub user_bound_scale: HighsInt,
    /// HiGHS: HighsOptions.h `blend_multi_objectives` (record lines 1424-1427: default true)
    pub blend_multi_objectives: bool,
    /// HiGHS: HighsOptions.h `use_warm_start` (record lines 1476-1479: default true)
    pub use_warm_start: bool,
    pub log_file: String,

    // IPM / IPX (HighsOptions.h, IPX; read by ipm/ipx_wrapper.rs)
    pub ipm_iteration_limit: HighsInt,
    pub ipm_optimality_tolerance: f64,
    pub start_crossover_tolerance: f64,
    pub ipx_dualize_strategy: HighsInt,
    pub run_centring: bool,
    pub max_centring_steps: HighsInt,
    pub centring_ratio_tolerance: f64,

    // PDLP (HighsOptions.h 403-410, 1339-1385). `pdlp_features_off`
    // has no OptionRecord upstream: it stays at its struct default 0.
    pub pdlp_features_off: HighsInt,
    pub pdlp_iteration_limit: HighsInt,
    pub pdlp_scaling_mode: HighsInt,
    pub pdlp_ruiz_iterations: HighsInt,
    pub pdlp_restart_strategy: HighsInt,
    pub pdlp_cupdlpc_restart_method: HighsInt,
    pub pdlp_step_size_strategy: HighsInt,
    pub pdlp_optimality_tolerance: f64,

    // HiPO (HighsOptions.h 397-400, 1312-1336). The struct
    // initialiser has `hipo_block_size(0)`; the OptionRecord default (128)
    // is what a Highs instance sees, so 128 is used here like the other
    // registered options.
    pub hipo_system: String,
    /// rusolve: what an interior point failure does: "off" is HiGHS
    /// (HiPO then IPX restart, then simplex clean-up); not a HiGHS option
    pub lp_fold: bool,
    pub lp_crossover_ipm_crash: bool,
    pub lp_crossover_batch_push: bool,
    /// rusolve: dual simplex from a non-logical basis uses Devex pricing
    /// instead of computing exact dual steepest edge weights (one BTRAN per row)
    pub lp_dual_devex_start: bool,
    /// rusolve: after postsolve, a nonbasic variable whose recovered value
    /// lies at the other bound than its status is switched to that bound
    pub lp_postsolve_snap_basis: bool,
    /// rusolve: HiPO support-reduction phase before crossover
    pub lp_ipm_support_reduction: bool,
    /// rusolve: race dual simplex against HiPO + crossover when threads >= 2
    pub lp_concurrent: bool,
    /// rusolve: pick ipm over simplex from structure, 1 thread
    pub lp_solver_auto: bool,
    /// rusolve: the clean-up solve after postsolve treats columns whose
    /// range is within the primal feasibility tolerance as fixed at their lower bound
    pub lp_postsolve_fix_tiny: bool,
    pub lp_fold_basis: bool,
    pub lp_crossover_ipx_refine: bool,
    pub lp_ipx_no_drop: bool,
    pub lp_ipx_soft_drop: bool,
    pub presolve_dep_eq_fixed_budget: bool,
    pub lp_dualize_auto: bool,
    /// rusolve: packed dgemm and grouped dtrsm in the HiPO factorisation.
    pub ipm_hipo_fast_blas: bool,
    /// rusolve: replay the recorded permutation of the HiPO matrix structure instead of redoing it every factorisation (bit-identical).
    pub ipm_hipo_cache_perm: bool,
    /// rusolve: build the HiPO normal-equations matrix on several threads, row blocks with private scatter arrays (bit-identical).
    pub ipm_hipo_par_ne: bool,
    /// rusolve: stop HiPO when dinf is flat while mu collapses and hand the LP to IPX.
    pub ipm_hipo_stall_detect: bool,
    /// rusolve: AMD first, Metis only above a flops threshold.
    pub ipm_hipo_amd_first: bool,
    pub ipm_hipo_dense_cols: bool,
    /// rusolve: dense columns of A out of the HiPO normal-equations factor, recovered by a Woodbury (Schur complement) update as CG preconditioner.
    pub ipm_hipo_dense_schur: bool,
    pub ipm_hipo_dualize: bool,
    pub lp_ipm_failure_policy: bool,
    pub hipo_parallel_type: String,
    pub hipo_ordering: String,
    pub hipo_block_size: HighsInt,
    /// Non-upstream: models `HighsExternalApi::isAvailable<HighsExtras::hipo>()`
    /// (HiGHS built with the HiPO extras). `false` keeps every dispatch site
    /// on the non-HiPO path; the HiPO tests set it. Not a registered option.
    pub hipo_available: bool,

    // QP / IIS / iCrash (HighsOptions.h)
    pub qp_allow_hot_start: bool,
    pub qp_iteration_limit: HighsInt,
    pub qp_nullspace_limit: HighsInt,
    pub qp_regularization_value: f64,
    pub iis_strategy: HighsInt,
    pub iis_time_limit: f64,
    pub icrash: bool,
    pub icrash_dualize: bool,
    pub icrash_strategy: String,
    pub icrash_starting_weight: f64,
    pub icrash_iterations: HighsInt,
    pub icrash_approx_iter: HighsInt,
    pub icrash_exact: bool,
    pub icrash_breakpoints: bool,

    // Tolerances
    pub primal_feasibility_tolerance: f64,
    pub dual_feasibility_tolerance: f64,
    pub kkt_tolerance: f64,
    pub primal_residual_tolerance: f64,
    pub dual_residual_tolerance: f64,
    pub optimality_tolerance: f64,
    pub infinite_cost: f64,
    pub infinite_bound: f64,
    pub small_matrix_value: f64,
    pub large_matrix_value: f64,

    // Simplex
    pub simplex_strategy: HighsInt,
    pub simplex_scale_strategy: HighsInt,
    pub simplex_dual_edge_weight_strategy: HighsInt,
    pub simplex_primal_edge_weight_strategy: HighsInt,
    pub simplex_iteration_limit: HighsInt,
    pub simplex_update_limit: HighsInt,
    pub simplex_min_concurrency: HighsInt,
    pub simplex_max_concurrency: HighsInt,
    pub simplex_price_strategy: HighsInt,
    pub simplex_unscaled_solution_strategy: HighsInt,
    pub dual_simplex_cost_perturbation_multiplier: f64,
    pub primal_simplex_bound_perturbation_multiplier: f64,
    pub dual_simplex_pivot_growth_tolerance: f64,
    pub dual_steepest_edge_weight_error_tolerance: f64,
    pub dual_steepest_edge_weight_log_error_threshold: f64,
    pub rebuild_refactor_solution_error_tolerance: f64,
    pub max_dual_simplex_cleanup_level: HighsInt,
    pub max_dual_simplex_phase1_cleanup_level: HighsInt,
    pub no_unnecessary_rebuild_refactor: bool,
    pub less_infeasible_DSE_check: bool,
    pub less_infeasible_DSE_choose_row: bool,
    pub allowed_cost_scale_factor: HighsInt,
    pub cost_scale_factor: HighsInt,
    pub allowed_matrix_scale_factor: HighsInt,
    pub simplex_dualize_strategy: HighsInt,
    pub simplex_permute_strategy: HighsInt,

    // Factorization
    pub factor_pivot_threshold: f64,
    pub factor_pivot_tolerance: f64,

    // Logging
    pub output_flag: bool,
    pub log_to_console: bool,
    pub log_dev_level: HighsInt,
    pub timeless_log: bool,
    /// HiGHS: HighsOptions::log_options. In C++ its flags point at
    /// `output_flag`/`log_to_console`/`log_dev_level`; here they are copies.
    /// After changing one of those three fields call [`HighsOptions::sync_log_options`].
    pub log_options: HighsLogOptions,
}

/// `h16_*` options [`HighsOptions::recommended`] leaves off. Empty, measured on our
/// benchmark set: the default is the HiGHS 1.16.0 behaviour. h16_zero_cost_singleton and h16_watch_link were here as
/// measured losses (pg 4.49 s against 38.4 s with zero_cost_singleton), but both are
/// genuine 1.16.0 behaviour and the faithful profile is
/// 0.983x time [0.942, 1.019], 129 vs 127 solved, against the profile without them.
pub const H16_RECOMMENDED_OFF: [&str; 0] = [];

/// Non-`h16_*` MIP options [`HighsOptions::recommended`] turns on: every `mip_*` /
/// `mip_*` 1.16.0 PR (each of those PRs is an ancestor of v1.16.0), so the default
/// profile is HiGHS 1.16.0. History of the set (measured on our benchmark set):
///
/// - #3213, #3244 first.
/// - The 1.16.0 PRs that do not change the search: together 0.995x time
///   [0.944, 1.034]; neos-787933 0.37x from #3404 alone.
/// - The clique trio (#3286/#3318/#3332), the #3312 presolve clique table and
///   Fourier-Motzkin presolve (#3218/#3250/#3255): 0.953x [0.888, 1.015], 213 vs 212
///   solved.
/// - #3264, #2990, #3209 (aggregate implications, machine-scheduling separator, keep
///   implications) and #3535 (64-bit symmetry hashes), with the two `h16_*` options above and
///   without the SCIP conflict-resolution and dynamic-symmetry options: the faithful 1.16.0
///   profile, 0.983x [0.942, 1.019], 129 vs 127 solved; uccase12 0.59x, neos-827175 0.74x,
///   krka 0.75x. With those SCIP options on top it is 0.982x [0.930, 1.029], so
///   the SCIP options are not default: they move the search off 1.16.0 for no measured
///   gain.
///
/// presolve_fm_obj_reformulation is not in the set: HiGHS 1.16.0 defaults it to false, and
/// turia runs 33.5 s with it against 13.4 s without.
pub const RECOMMENDED_MIP_ON: [&str; 25] = [
    "mip_sb_budget_pre_restart",    // HiGHS 1.16.0 #3213
    "mip_dual_fix_probing",         // HiGHS 1.16.0 #3244
    "mip_cdouble_act_contrib",      // HiGHS 1.16.0 #3184
    "mip_dfs_node_limit",           // HiGHS 1.16.0 #3445
    "mip_dual_subst_limits",        // HiGHS 1.16.0 #3404
    "mip_keep_branch_reason",       // HiGHS 1.16.0 #3357
    "mip_lurking_bounds_robust",    // HiGHS 1.16.0 #3481, #3474
    "mip_path_same_bds",            // HiGHS 1.16.0 #3179
    "mip_ranged_dual_implied",      // HiGHS 1.16.0 #3437
    "mip_reduced_stack_fast",       // HiGHS 1.16.0 #3336 (b)
    "mip_redundant_vbd",            // HiGHS 1.16.0 #3177
    "mip_strengthen_small_contrib", // HiGHS 1.16.0 #3468
    "mip_cdouble_proofs",           // HiGHS 1.16.0 #3223
    "mip_fix_last_active",          // HiGHS 1.16.0 #3220
    "mip_clique_normalise",         // HiGHS 1.16.0 #3286
    "mip_coef_change_impl_bounds",  // HiGHS 1.16.0 #3318
    "mip_dual_side_lost",           // HiGHS 1.16.0 #3332
    "mip_presolve_clique_table",    // HiGHS 1.16.0 #3312
    "mip_presolve_row_append",      // HiGHS 1.16.0 #3218, #3250
    "mip_fm_presolve",              // HiGHS 1.16.0 #3255
    "mip_aggregate_implications",   // HiGHS 1.16.0 #3264
    "mip_machine_sched_sepa",       // HiGHS 1.16.0 #2990
    "mip_keep_implications",        // HiGHS 1.16.0 #3209
    "mip_sym_hash64",               // HiGHS 1.16.0 #3535
    // rusolve, no HiGHS PR: exits the 1.16.0 strengthenInequalities cover loop where it
    // would spin forever; search-identical
    // with it off on all 44 instances of our benchmark set, so it only changes runs where 1.16.0 hangs.
    "mip_strengthen_cover_guard",
];

impl HighsOptions {
    /// The best-measured profile, the CLI default: [`HighsOptions::default`]
    /// (HiGHS 1.15.1 semantics, every `h16_*` option off) plus every `h16_*` option
    /// except [`H16_RECOMMENDED_OFF`], plus [`RECOMMENDED_MIP_ON`]: HiGHS 1.16.0.
    /// Sub-solvers built internally with `Highs::default()` (C++ `Highs sub;`) get
    /// the version flags through h16_sub_solver_flags.
    pub fn recommended() -> Self {
        let mut options = Self::default();
        options.set_all_h16(true);
        options.set_recommended_mip(true);
        options
    }

    /// Sets every option in [`RECOMMENDED_MIP_ON`] to `on`.
    pub fn set_recommended_mip(&mut self, on: bool) {
        self.mip_sb_budget_pre_restart = on;
        self.mip_dual_fix_probing = on;
        self.mip_cdouble_act_contrib = on;
        self.mip_dfs_node_limit = on;
        self.mip_dual_subst_limits = on;
        self.mip_keep_branch_reason = on;
        self.mip_lurking_bounds_robust = on;
        self.mip_path_same_bds = on;
        self.mip_ranged_dual_implied = on;
        self.mip_reduced_stack_fast = on;
        self.mip_redundant_vbd = on;
        self.mip_strengthen_small_contrib = on;
        self.mip_cdouble_proofs = on;
        self.mip_fix_last_active = on;
        self.mip_clique_normalise = on;
        self.mip_coef_change_impl_bounds = on;
        self.mip_dual_side_lost = on;
        self.mip_presolve_clique_table = on;
        self.mip_presolve_row_append = on;
        self.mip_fm_presolve = on;
        self.mip_aggregate_implications = on;
        self.mip_machine_sched_sepa = on;
        self.mip_keep_implications = on;
        self.mip_sym_hash64 = on;
        self.mip_strengthen_cover_guard = on;
    }

    /// Copies every `h16_*`, `mip_*` and `mip_*` flag (the 1.16.0 code switches) from `src`. Used by
    /// the MIP's internal solvers under `h16_sub_solver_flags`: HiGHS 1.16.0 runs those
    /// solvers on 1.16.0 code although their options start at the defaults.
    pub fn copy_version_flags_from(&mut self, src: &HighsOptions) {
        self.h16_ac_centring = src.h16_ac_centring;
        self.h16_ac_fix = src.h16_ac_fix;
        self.h16_add_no_modify_on_error = src.h16_add_no_modify_on_error;
        self.h16_agg_term_zero = src.h16_agg_term_zero;
        self.h16_assess_small = src.h16_assess_small;
        self.h16_big_act = src.h16_big_act;
        self.h16_bound_chain = src.h16_bound_chain;
        self.h16_bound_range_inf = src.h16_bound_range_inf;
        self.h16_cdouble = src.h16_cdouble;
        self.h16_clique_links = src.h16_clique_links;
        self.h16_col_stuffing = src.h16_col_stuffing;
        self.h16_cost_perturb = src.h16_cost_perturb;
        self.h16_crossover_x = src.h16_crossover_x;
        self.h16_cupdlp_time = src.h16_cupdlp_time;
        self.h16_cutpool_age = src.h16_cutpool_age;
        self.h16_dep_eq_time_cap = src.h16_dep_eq_time_cap;
        self.h16_dominated_cols = src.h16_dominated_cols;
        self.h16_dual_free_list = src.h16_dual_free_list;
        self.h16_dual_postsolve_rules = src.h16_dual_postsolve_rules;
        self.h16_dup_row_undo = src.h16_dup_row_undo;
        self.h16_dupcol_warn = src.h16_dupcol_warn;
        self.h16_enum_sig64 = src.h16_enum_sig64;
        self.h16_eq_row_adds = src.h16_eq_row_adds;
        self.h16_feasrelax_sense = src.h16_feasrelax_sense;
        self.h16_fj = src.h16_fj;
        self.h16_frontier_merge = src.h16_frontier_merge;
        self.h16_hash_common = src.h16_hash_common;
        self.h16_hfactor_pivot = src.h16_hfactor_pivot;
        self.h16_hipdlp_max = src.h16_hipdlp_max;
        self.h16_hipo_display = src.h16_hipo_display;
        self.h16_hipo_empty_rows = src.h16_hipo_empty_rows;
        self.h16_impl_eq = src.h16_impl_eq;
        self.h16_implint_bounds = src.h16_implint_bounds;
        self.h16_int_bound_round = src.h16_int_bound_round;
        self.h16_int_bound_rounding = src.h16_int_bound_rounding;
        self.h16_int_round_assert_tol = src.h16_int_round_assert_tol;
        self.h16_ipm_gap_offset = src.h16_ipm_gap_offset;
        self.h16_ipx_free_rows = src.h16_ipx_free_rows;
        self.h16_ipx_scaled = src.h16_ipx_scaled;
        self.h16_ipx_stall = src.h16_ipx_stall;
        self.h16_kkt_effective_costs = src.h16_kkt_effective_costs;
        self.h16_max_duals = src.h16_max_duals;
        self.h16_mip_start = src.h16_mip_start;
        self.h16_mip_workers = src.h16_mip_workers;
        self.h16_ninfmin = src.h16_ninfmin;
        self.h16_nla_debug_gate = src.h16_nla_debug_gate;
        self.h16_obj_cdouble = src.h16_obj_cdouble;
        self.h16_parallel_rows = src.h16_parallel_rows;
        self.h16_postsolve_basis_check = src.h16_postsolve_basis_check;
        self.h16_presolve_rules = src.h16_presolve_rules;
        self.h16_row_duals = src.h16_row_duals;
        self.h16_fourier_motzkin = src.h16_fourier_motzkin;
        self.h16_qp_bounded_start = src.h16_qp_bounded_start;
        self.h16_qp_convexity = src.h16_qp_convexity;
        self.h16_ray_simplex = src.h16_ray_simplex;
        self.h16_rebuild_reason = src.h16_rebuild_reason;
        self.h16_repair_lp = src.h16_repair_lp;
        self.h16_row_presolve = src.h16_row_presolve;
        self.h16_row_redundant = src.h16_row_redundant;
        self.h16_scaling_report = src.h16_scaling_report;
        self.h16_singleton_loops = src.h16_singleton_loops;
        self.h16_subst_int = src.h16_subst_int;
        self.h16_sym_refine = src.h16_sym_refine;
        self.h16_unresolved_unbounded_or_infeasible = src.h16_unresolved_unbounded_or_infeasible;
        self.h16_user_scale_offset = src.h16_user_scale_offset;
        self.h16_user_sol_obj = src.h16_user_sol_obj;
        self.h16_user_solution_assess = src.h16_user_solution_assess;
        self.h16_watch_link = src.h16_watch_link;
        self.h16_zero_cost_singleton = src.h16_zero_cost_singleton;
        self.h16_sub_solver_flags = src.h16_sub_solver_flags;
        self.mip_sb_budget_pre_restart = src.mip_sb_budget_pre_restart;
        self.mip_dfs_node_limit = src.mip_dfs_node_limit;
        self.mip_reduced_stack_fast = src.mip_reduced_stack_fast;
        self.mip_cdouble_act_contrib = src.mip_cdouble_act_contrib;
        self.mip_keep_branch_reason = src.mip_keep_branch_reason;
        self.mip_redundant_vbd = src.mip_redundant_vbd;
        self.mip_path_same_bds = src.mip_path_same_bds;
        self.mip_strengthen_small_contrib = src.mip_strengthen_small_contrib;
        self.mip_lurking_bounds_robust = src.mip_lurking_bounds_robust;
        self.mip_dual_subst_limits = src.mip_dual_subst_limits;
        self.mip_ranged_dual_implied = src.mip_ranged_dual_implied;
        self.mip_machine_sched_sepa = src.mip_machine_sched_sepa;
        self.mip_dual_fix_probing = src.mip_dual_fix_probing;
        self.mip_aggregate_implications = src.mip_aggregate_implications;
        self.mip_clique_normalise = src.mip_clique_normalise;
        self.mip_coef_change_impl_bounds = src.mip_coef_change_impl_bounds;
        self.mip_dual_side_lost = src.mip_dual_side_lost;
        self.mip_cdouble_proofs = src.mip_cdouble_proofs;
        self.mip_fix_last_active = src.mip_fix_last_active;
        self.mip_keep_implications = src.mip_keep_implications;
        self.mip_presolve_clique_table = src.mip_presolve_clique_table;
        self.mip_presolve_row_append = src.mip_presolve_row_append;
        self.mip_fm_presolve = src.mip_fm_presolve;
        self.mip_sym_hash64 = src.mip_sym_hash64;
    }

    /// Sets every `h16_*` option to `on`.
    pub fn set_all_h16(&mut self, on: bool) {
        self.h16_ac_centring = on;
        self.h16_ac_fix = on;
        self.h16_add_no_modify_on_error = on;
        self.h16_agg_term_zero = on;
        self.h16_assess_small = on;
        self.h16_big_act = on;
        self.h16_bound_chain = on;
        self.h16_bound_range_inf = on;
        self.h16_cdouble = on;
        self.h16_clique_links = on;
        self.h16_col_stuffing = on;
        self.h16_cost_perturb = on;
        self.h16_crossover_x = on;
        self.h16_cupdlp_time = on;
        self.h16_cutpool_age = on;
        self.h16_dep_eq_time_cap = on;
        self.h16_dominated_cols = on;
        self.h16_dual_free_list = on;
        self.h16_dual_postsolve_rules = on;
        self.h16_dup_row_undo = on;
        self.h16_dupcol_warn = on;
        self.h16_enum_sig64 = on;
        self.h16_eq_row_adds = on;
        self.h16_feasrelax_sense = on;
        self.h16_fj = on;
        self.h16_frontier_merge = on;
        self.h16_hash_common = on;
        self.h16_hfactor_pivot = on;
        self.h16_hipdlp_max = on;
        self.h16_hipo_display = on;
        self.h16_hipo_empty_rows = on;
        self.h16_impl_eq = on;
        self.h16_implint_bounds = on;
        self.h16_int_bound_round = on;
        self.h16_int_bound_rounding = on;
        self.h16_int_round_assert_tol = on;
        self.h16_ipm_gap_offset = on;
        self.h16_ipx_free_rows = on;
        self.h16_ipx_scaled = on;
        self.h16_ipx_stall = on;
        self.h16_kkt_effective_costs = on;
        self.h16_max_duals = on;
        self.h16_mip_start = on;
        self.h16_mip_workers = on;
        self.h16_ninfmin = on;
        self.h16_nla_debug_gate = on;
        self.h16_obj_cdouble = on;
        self.h16_parallel_rows = on;
        self.h16_postsolve_basis_check = on;
        self.h16_presolve_rules = on;
        self.h16_row_duals = on;
        self.h16_fourier_motzkin = on;
        self.h16_qp_bounded_start = on;
        self.h16_qp_convexity = on;
        self.h16_ray_simplex = on;
        self.h16_rebuild_reason = on;
        self.h16_repair_lp = on;
        self.h16_row_presolve = on;
        self.h16_row_redundant = on;
        self.h16_scaling_report = on;
        self.h16_singleton_loops = on;
        self.h16_subst_int = on;
        self.h16_sym_refine = on;
        self.h16_unresolved_unbounded_or_infeasible = on;
        self.h16_user_scale_offset = on;
        self.h16_user_sol_obj = on;
        self.h16_user_solution_assess = on;
        self.h16_watch_link = on;
        self.h16_zero_cost_singleton = on;
        self.h16_sub_solver_flags = on;
    }
}

impl Default for HighsOptions {
    fn default() -> Self {
        Self {
            time_limit: kHighsInf,
            objective_bound: kHighsInf,
            objective_target: -kHighsInf,
            threads: 0,
            random_seed: 0,
            parallel: kHighsChooseString.to_string(),
            highs_debug_level: kHighsDebugLevelNone,
            highs_analysis_level: kHighsAnalysisLevelNone,
            allow_unbounded_or_infeasible: false,

            presolve: kHighsChooseString.to_string(),
            solver: kHighsChooseString.to_string(),
            run_crossover: kHighsOnString.to_string(),
            use_implied_bounds_from_presolve: false,
            lp_presolve_requires_basis_postsolve: true,
            presolve_reduction_limit: -1,
            restart_presolve_reduction_limit: -1,
            presolve_substitution_maxfillin: 10,
            presolve_rule_off: 0,
            presolve_rule_test: 0,
            presolve_rule_logging: false,
            presolve_remove_slacks: false,
            presolve_pivot_threshold: 0.01,
            mip_feasibility_tolerance: kDefaultMipTolerance,

            mip_detect_symmetry: true,
            mip_allow_restart: true,
            mip_max_nodes: kHighsIInf,
            mip_max_stall_nodes: kHighsIInf,
            mip_max_start_nodes: 500,
            mip_debug_solution_file: String::new(),
            mip_improving_solution_save: false,
            mip_improving_solution_report_sparse: false,
            mip_improving_solution_file: String::new(),
            mip_root_presolve_only: false,
            mip_lifting_for_probing: -1,
            mip_max_leaves: kHighsIInf,
            mip_max_improving_sols: kHighsIInf,
            mip_lp_age_limit: 10,
            mip_pool_age_limit: 30,
            mip_pool_soft_limit: 10000,
            mip_pscost_minreliable: 8,
            mip_min_cliquetable_entries_for_parallelism: 100000,
            mip_report_level: 1,
            mip_heuristic_effort: 0.05,
            mip_heuristic_run_feasibility_jump: true,
            mip_heuristic_run_rins: true,
            mip_heuristic_run_rens: true,
            mip_heuristic_run_root_reduced_cost: true,
            mip_heuristic_run_zi_round: false,
            mip_heuristic_run_shifting: false,
            mip_allow_cut_separation_at_nodes: true,
            mip_rel_gap: 1e-4,
            mip_abs_gap: 1e-6,
            mip_min_logging_interval: 5.0,
            mip_lp_solver: kHighsChooseString.to_string(),
            mip_ipm_solver: kHighsChooseString.to_string(),
            mip_search_simulate_concurrency: false,
            mip_sb_budget_pre_restart: false,
            mip_dfs_node_limit: false,
            mip_reduced_stack_fast: false,
            mip_cdouble_act_contrib: false,
            mip_keep_branch_reason: false,
            mip_redundant_vbd: false,
            mip_path_same_bds: false,
            mip_strengthen_small_contrib: false,
            mip_lurking_bounds_robust: false,
            mip_dual_subst_limits: false,
            mip_ranged_dual_implied: false,
            mip_machine_sched_sepa: false,
            mip_dual_fix_probing: false,
            mip_aggregate_implications: false,
            mip_clique_normalise: false,
            mip_coef_change_impl_bounds: false,
            mip_dual_side_lost: false,
            mip_cdouble_proofs: false,
            mip_fix_last_active: false,
            mip_keep_implications: false,
            mip_presolve_clique_table: false,
            mip_presolve_row_append: false,
            mip_fm_presolve: false,
            mip_sym_hash64: false,
            // rusolve backport: HiGHS 1.16.0 options, default false
            h16_ac_centring: false,
            h16_ac_fix: false,
            h16_add_no_modify_on_error: false,
            h16_agg_term_zero: false,
            h16_assess_small: false,
            h16_big_act: false,
            h16_bound_chain: false,
            h16_bound_range_inf: false,
            h16_cdouble: false,
            h16_clique_links: false,
            h16_col_stuffing: false,
            h16_cost_perturb: false,
            h16_crossover_x: false,
            h16_cupdlp_time: false,
            h16_cutpool_age: false,
            h16_dep_eq_time_cap: false,
            h16_dominated_cols: false,
            h16_dual_free_list: false,
            h16_dual_postsolve_rules: false,
            h16_dup_row_undo: false,
            h16_dupcol_warn: false,
            h16_enum_sig64: false,
            h16_eq_row_adds: false,
            h16_feasrelax_sense: false,
            h16_fj: false,
            h16_frontier_merge: false,
            h16_hash_common: false,
            h16_hfactor_pivot: false,
            h16_hipdlp_max: false,
            h16_hipo_display: false,
            h16_hipo_empty_rows: false,
            h16_impl_eq: false,
            h16_implint_bounds: false,
            h16_int_bound_round: false,
            h16_int_bound_rounding: false,
            h16_int_round_assert_tol: false,
            h16_ipm_gap_offset: false,
            h16_ipx_free_rows: false,
            h16_ipx_scaled: false,
            h16_ipx_stall: false,
            h16_kkt_effective_costs: false,
            h16_max_duals: false,
            h16_mip_start: false,
            h16_mip_workers: false,
            h16_ninfmin: false,
            h16_nla_debug_gate: false,
            h16_obj_cdouble: false,
            h16_parallel_rows: false,
            h16_postsolve_basis_check: false,
            h16_presolve_rules: false,
            presolve_mode: kPresolveModeDefault,
            presolve_fm_obj_reformulation: false,
            h16_row_duals: false,
            h16_fourier_motzkin: false,
            h16_qp_bounded_start: false,
            h16_qp_convexity: false,
            h16_ray_simplex: false,
            h16_rebuild_reason: false,
            h16_repair_lp: false,
            h16_row_presolve: false,
            h16_row_redundant: false,
            h16_scaling_report: false,
            h16_singleton_loops: false,
            h16_subst_int: false,
            h16_sym_refine: false,
            h16_unresolved_unbounded_or_infeasible: false,
            h16_user_scale_offset: false,
            h16_user_sol_obj: false,
            h16_user_solution_assess: false,
            h16_watch_link: false,
            h16_zero_cost_singleton: false,
            h16_sub_solver_flags: false,
            mip_conflict_resolution: false,
            mip_conflict_resolution_tree_only: false,
            mip_symmetry_dynamic: false,
            mip_cut_selection_hybrid: false,
            mip_knapsack_cover_sepa: false,
            mip_dff_separator: false,
            mip_tableau_loosen: false,
            mip_tableau_loosen_auto: false,
            mip_ti_precedence_sepa: true,
            mip_ti_list_sched: true,
            mip_skip_analytic_centre: false,
            mip_netdesign_sepa: true,
            mip_netdesign_heur: true,
            mip_netdesign_cutset: false,
            mip_cutoff_bigm: false,
            mip_node_bigm: false,
            mip_node_bigm_rows: false,
            mip_node_bigm_min_ratio: 0.001,
            mip_node_bigm_active: 0.0,
            mip_cont_prop_improve: 0.05,
            mip_pscost_variance_tol: 0.0,
            mip_rens_min_fixing_rate: 0.0,
            mip_heuristic_run_oneopt: false,
            mip_heuristic_run_crossover: false,
            mip_heuristic_run_fpr: false,
            mip_heuristic_run_local_mip: false,
            mip_strengthen_cover_guard: false,
            mip_ls_gub_swap: false,
            mip_ls_int_lp_repair: false,
            mip_alns: false,
            mip_presolve_eq_subst: false,
            mip_presolve_clique_varbound: false,
            mip_presolve_xor_gf2: true,
            mip_orbitope_root_fix: true,
            mip_orbitope_full_fix: false,
            mip_shift_and_propagate: false,
            mip_proximity: false,
            mip_improve_early: false,
            mip_primal_stall: false,
            mip_improve_early_work: 20.0,
            mip_guided_dive: false,
            mip_local_branching: false,
            mip_cdcl_dive: false,
            mip_root_bound_stall: true,
            mip_cdcl_uip: false,
            mip_zero_obj_cdcl: true,
            mip_decompose: false,
            mip_lattice: true,
            mip_activation_dff: false,
            mip_fj_native: false,
            mip_first_inc: false,
            mip_first_inc_root: false,
            mip_lattice_kernel: true,
            mip_presolve_probing_cap: false,
            mip_presolve_domcol_full: false,
            mip_presolve_probing_useless: 0,
            mip_conflict_off_mask: 0,
            mip_redcost_proof_uncapped: false,
            mip_root_stall_scip: false,
            mip_heur_lp_quota: true,
            mip_racing_threads: 0,
            mip_root_lp_concurrent: false,
            solve_relaxation: false,
            user_objective_scale: 0,
            user_bound_scale: 0,
            blend_multi_objectives: true,
            use_warm_start: true,
            log_file: String::new(),

            ipm_iteration_limit: kHighsIInf,
            ipm_optimality_tolerance: kDefaultIpmTolerance,
            start_crossover_tolerance: kDefaultIpmTolerance,
            ipx_dualize_strategy: kIpxDualizeStrategyLukas,
            run_centring: false,
            max_centring_steps: 5,
            centring_ratio_tolerance: 100.0,

            pdlp_features_off: 0,
            pdlp_iteration_limit: kHighsIInf,
            pdlp_scaling_mode: kPdlpScalingRuiz + kPdlpScalingPC,
            pdlp_ruiz_iterations: 10,
            pdlp_restart_strategy: kPdlpRestartStrategyAdaptive,
            pdlp_cupdlpc_restart_method: 1,
            pdlp_step_size_strategy: kPdlpStepSizeStrategyAdaptive,
            pdlp_optimality_tolerance: kDefaultKktTolerance,

            hipo_system: kHighsChooseString.to_string(),
            lp_fold: true,
            lp_crossover_ipm_crash: false,
            lp_crossover_batch_push: false,
            lp_dual_devex_start: false,
            lp_postsolve_snap_basis: false,
            lp_ipm_support_reduction: false,
            lp_concurrent: true,
            lp_solver_auto: true,
            lp_postsolve_fix_tiny: true,
            lp_fold_basis: false,
            lp_crossover_ipx_refine: true,
            lp_ipx_no_drop: false,
            lp_ipx_soft_drop: true,
            presolve_dep_eq_fixed_budget: false,
            lp_dualize_auto: true,
            ipm_hipo_fast_blas: true,
            ipm_hipo_cache_perm: true,
            ipm_hipo_par_ne: true,
            ipm_hipo_stall_detect: false,
            ipm_hipo_amd_first: true,
            ipm_hipo_dense_cols: false,
            ipm_hipo_dense_schur: false,
            ipm_hipo_dualize: true,
            lp_ipm_failure_policy: true,
            hipo_parallel_type: kHipoBothString.to_string(),
            hipo_ordering: kHighsChooseString.to_string(),
            hipo_block_size: 128,
            hipo_available: false,

            qp_allow_hot_start: false,
            qp_iteration_limit: kHighsIInf,
            qp_nullspace_limit: 4000,
            qp_regularization_value: kHessianRegularizationValue,
            iis_strategy: kIisStrategyDefault,
            iis_time_limit: kHighsInf,
            icrash: false,
            icrash_dualize: false,
            icrash_strategy: "ICA".to_string(),
            icrash_starting_weight: 1e-3,
            icrash_iterations: 30,
            icrash_approx_iter: 50,
            icrash_exact: false,
            icrash_breakpoints: false,

            primal_feasibility_tolerance: 1e-7,
            dual_feasibility_tolerance: 1e-7,
            kkt_tolerance: kDefaultKktTolerance,
            primal_residual_tolerance: kDefaultKktTolerance,
            dual_residual_tolerance: kDefaultKktTolerance,
            optimality_tolerance: kDefaultKktTolerance,
            infinite_cost: 1e20,
            infinite_bound: 1e20,
            small_matrix_value: 1e-9,
            large_matrix_value: 1e15,

            simplex_strategy: kSimplexStrategyDual,
            simplex_scale_strategy: kSimplexScaleStrategyEquilibration,
            simplex_dual_edge_weight_strategy: kSimplexEdgeWeightStrategyChoose,
            simplex_primal_edge_weight_strategy: kSimplexEdgeWeightStrategyChoose,
            simplex_iteration_limit: kHighsIInf,
            simplex_update_limit: 5000,
            simplex_min_concurrency: 1,
            simplex_max_concurrency: kSimplexConcurrencyLimit,
            simplex_price_strategy: kSimplexPriceStrategyRowSwitchColSwitch,
            simplex_unscaled_solution_strategy: kSimplexUnscaledSolutionStrategyRefine,
            dual_simplex_cost_perturbation_multiplier: 1.0,
            primal_simplex_bound_perturbation_multiplier: 1.0,
            dual_simplex_pivot_growth_tolerance: 1e-9,
            dual_steepest_edge_weight_error_tolerance: kHighsInf,
            dual_steepest_edge_weight_log_error_threshold: 1e1,
            rebuild_refactor_solution_error_tolerance: 1e-8,
            max_dual_simplex_cleanup_level: 1,
            max_dual_simplex_phase1_cleanup_level: 2,
            no_unnecessary_rebuild_refactor: true,
            less_infeasible_DSE_check: true,
            less_infeasible_DSE_choose_row: true,
            allowed_cost_scale_factor: 0,
            cost_scale_factor: 0,
            allowed_matrix_scale_factor: kDefaultAllowedMatrixPow2Scale,
            simplex_dualize_strategy: kHighsOptionOff,
            simplex_permute_strategy: kHighsOptionOff,

            factor_pivot_threshold: 0.1,
            factor_pivot_tolerance: 1e-10,

            output_flag: true,
            log_to_console: true,
            log_dev_level: kHighsLogDevLevelNone,
            timeless_log: false,
            log_options: HighsLogOptions::default(),
        }
    }
}

impl HighsOptions {
    /// Copy `output_flag`, `log_to_console`, `log_dev_level` into
    /// `log_options` (C++ shares them by pointer, so they are always in sync).
    pub fn sync_log_options(&mut self) {
        self.log_options.output_flag = self.output_flag;
        self.log_options.log_to_console = self.log_to_console;
        self.log_options.log_dev_level = self.log_dev_level;
    }
}

use crate::io::highs_io::kHighsLogDevLevelNone;

/// HiGHS: HighsOptions.h::kSimplexString
pub const kSimplexString: &str = "simplex";
/// HiGHS: HighsOptions.h::kIpmString
pub const kIpmString: &str = "ipm";
/// HiGHS: HighsOptions.h::kHipoString
pub const kHipoString: &str = "hipo";
/// HiGHS: HighsOptions.h::kIpxString
pub const kIpxString: &str = "ipx";
/// HiGHS: HighsOptions.h::kPdlpString
pub const kPdlpString: &str = "pdlp";
/// HiGHS: HighsOptions.h::kHiPdlpString
pub const kHiPdlpString: &str = "hipdlp";

/// HiGHS: HighsOptions.h::kSolverString
pub const kSolverString: &str = "solver";
/// HiGHS: HighsOptions.h::kMipLpSolverString
pub const kMipLpSolverString: &str = "mip_lp_solver";
/// HiGHS: HighsOptions.h::kMipIpmSolverString
pub const kMipIpmSolverString: &str = "mip_ipm_solver";
/// HiGHS: HighsOptions.h::kHipoSystemString
pub const kHipoSystemString: &str = "hipo_system";
/// HiGHS: HighsOptions.h::kHipoAugmentedString
pub const kHipoAugmentedString: &str = "augmented";
/// HiGHS: HighsOptions.h::kHipoNormalEqString
pub const kHipoNormalEqString: &str = "normaleq";
/// HiGHS: HighsOptions.h::kHipoParallelString
pub const kHipoParallelString: &str = "hipo_parallel_type";
/// HiGHS: HighsOptions.h::kHipoTreeString
pub const kHipoTreeString: &str = "tree";
/// HiGHS: HighsOptions.h::kHipoNodeString
pub const kHipoNodeString: &str = "node";
/// HiGHS: HighsOptions.h::kHipoBothString
pub const kHipoBothString: &str = "both";
/// HiGHS: HighsOptions.h::kHipoOrderingString
pub const kHipoOrderingString: &str = "hipo_ordering";
/// HiGHS: HighsOptions.h::kHipoMetisString
pub const kHipoMetisString: &str = "metis";
/// HiGHS: HighsOptions.h::kHipoAmdString
pub const kHipoAmdString: &str = "amd";
/// HiGHS: HighsOptions.h::kHipoRcmString
pub const kHipoRcmString: &str = "rcm";

/// HiGHS: HighsSolve.cpp::useIpm (lines 682-684, 3 lines)
pub fn use_ipm(solver: &str) -> bool {
    solver == kIpmString || solver == kHipoString || solver == kIpxString
}

/// HiGHS: HighsSolve.cpp::usePdlp (lines 686-688, 3 lines)
pub fn use_pdlp(solver: &str) -> bool {
    solver == kPdlpString || solver == kHiPdlpString
}

/// HiGHS: HighsSolve.cpp::useHipo (lines 691-720, 30 lines)
///
/// `HighsExternalApi::isAvailable<HighsExtras::hipo>()` is
/// `options.hipo_available` (field doc). `lp` and `logging` are unused in
/// 1.15.1 too ("Later decide ... based on LP properties").
pub fn use_hipo(
    options: &HighsOptions,
    specific_solver_option: &str,
    _lp: &crate::lp_data::highs_lp::HighsLp,
    _logging: bool,
) -> bool {
    debug_assert!(
        specific_solver_option == kSolverString
            || specific_solver_option == kMipLpSolverString
            || specific_solver_option == kMipIpmSolverString
    );
    let specific_solver_option_value: &str = if specific_solver_option == kSolverString {
        &options.solver
    } else if specific_solver_option == kMipLpSolverString {
        &options.mip_lp_solver
    } else {
        &options.mip_ipm_solver
    };
    // In the MIP solver there are situations where IPM must be used
    let force_ipm = specific_solver_option == kMipIpmSolverString;
    let mut use_hipo = false;
    if specific_solver_option_value == kIpxString {
        use_hipo = false;
    } else if specific_solver_option_value == kIpmString
        || specific_solver_option_value == kHipoString
        || force_ipm
    {
        use_hipo = options.hipo_available;
    }
    if options.run_centring {
        use_hipo = false;
    }
    use_hipo
}

/// HiGHS: HighsOptions.cpp::solverValidForLp (lines 1169-1174, 6 lines)
pub fn solver_valid_for_lp(solver: &str) -> bool {
    solver == kHighsChooseString
        || solver == kSimplexString
        || solver == kIpmString
        || solver == kIpxString
        || solver == kHipoString
        || solver == kPdlpString
        || solver == kHiPdlpString
}
