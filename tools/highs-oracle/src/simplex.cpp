// End-to-end simplex traces through the observer hooks added by
// ../highs-oracle-hooks.patch. Hooks only read solver state.
#include <cmath>

#include "lp_data/HighsOptions.h"
#include "oracle.h"
#include "simplex/HEkk.h"
#include "simplex/HEkkDual.h"
#include "simplex/HEkkPrimal.h"
#include "simplex/HighsOracleHook.h"

// Declared friend of HEkkPrimal and HSimplexNla by the patch.
struct HighsOracleAccess {
  static const HFactor& factor(const HEkk& ekk) { return ekk.simplex_nla_.factor_; }
  // `pivot` = false before the first iteration of this HEkkPrimal: the pivot
  // members (move_in ... numericalTrouble) are then uninitialised (the
  // constructor does not set them), so they are omitted.
  static Obj primalIter(const HEkkPrimal& p, bool pivot) {
    Obj o;
    o.i("solve_phase", p.solve_phase).i("rebuild_reason", p.rebuild_reason)
        .i("variable_in", p.variable_in);
    if (pivot) o.i("move_in", p.move_in);
    o.i("row_out", p.row_out);
    if (pivot)
      o.i("variable_out", p.variable_out)
          .i("move_out", p.move_out).f("theta_dual", p.theta_dual)
          .f("theta_primal", p.theta_primal).f("value_in", p.value_in)
          .f("alpha_col", p.alpha_col).f("alpha_row", p.alpha_row)
          .f("numericalTrouble", p.numericalTrouble);
    o.i("num_flip_since_rebuild", p.num_flip_since_rebuild)
        .i("edge_weight_mode", (int)p.edge_weight_mode)
        .h("h_edge_weight", fnvF(p.edge_weight_));
    if (p.variable_in >= 0 && p.variable_in < (HighsInt)p.edge_weight_.size())
      o.f("edge_weight_in", p.edge_weight_[p.variable_in]);
    return o;
  }
  static Obj primalState(const HEkkPrimal& p) {
    Obj o;
    o.i("solve_phase", p.solve_phase).i("rebuild_reason", p.rebuild_reason)
        .i("edge_weight_mode", (int)p.edge_weight_mode)
        .f("primal_feasibility_tolerance", p.primal_feasibility_tolerance)
        .f("dual_feasibility_tolerance", p.dual_feasibility_tolerance)
        .i("num_free_col", p.num_free_col)
        .i("num_flip_since_rebuild", p.num_flip_since_rebuild)
        .h("h_edge_weight", fnvF(p.edge_weight_));
    return o;
  }
  static const HEkk& ekk(const HEkkPrimal& p) { return p.ekk_instance_; }
  static HighsInt rebuildReason(const HEkkPrimal& p) { return p.rebuild_reason; }
};

namespace {

struct TraceState {
  Out* out = nullptr;
  HighsInt solve_call = -1;
  HighsInt dual_call = -1;
  HighsInt primal_call = -1;
  HighsInt last_iter = -1;
  HighsInt pending_reason = 0;  // rebuild_reason on entry to rebuild()
  bool primal_iterated = false;  // this HEkkPrimal::solve has iterated
  // snapshot mode
  EkkSnapshot* snap = nullptr;
};
TraceState g;

Obj optionsObj(const HighsOptions& options) {
  Obj o;
  for (const OptionRecord* r : options.records) {
    switch (r->type) {
      case HighsOptionType::kBool:
        o.b(r->name.c_str(), *static_cast<const OptionRecordBool*>(r)->value);
        break;
      case HighsOptionType::kInt:
        o.i(r->name.c_str(), *static_cast<const OptionRecordInt*>(r)->value);
        break;
      case HighsOptionType::kDouble:
        o.f(r->name.c_str(), *static_cast<const OptionRecordDouble*>(r)->value);
        break;
      case HighsOptionType::kString: {
        // Paths/files are environment, not solver state.
        if (r->name == "log_file" || r->name == "solution_file" ||
            r->name == "write_model_file" || r->name == "read_solution_file" ||
            r->name == "read_basis_file" || r->name == "write_basis_file" ||
            r->name == "write_presolved_model_file" || r->name == "write_iis_model_file" ||
            r->name == "glpsol_cost_row_location" || r->name == "solve_result_file")
          break;
        o.s(r->name.c_str(), *static_cast<const OptionRecordString*>(r)->value);
        break;
      }
    }
  }
  return o;
}

Obj lpObj(const HighsLp& lp) {
  Obj o;
  o.i("num_col", lp.num_col_).i("num_row", lp.num_row_).i("sense", (int)lp.sense_)
      .f("offset", lp.offset_).vf("col_cost", lp.col_cost_)
      .vf("col_lower", lp.col_lower_).vf("col_upper", lp.col_upper_)
      .vf("row_lower", lp.row_lower_).vf("row_upper", lp.row_upper_)
      .mat("a", lp.a_matrix_).b("is_scaled", lp.is_scaled_)
      .b("is_moved", lp.is_moved_);
  Obj s;
  s.i("strategy", lp.scale_.strategy).b("has_scaling", lp.scale_.has_scaling)
      .i("num_col", lp.scale_.num_col).i("num_row", lp.scale_.num_row)
      .f("cost", lp.scale_.cost).vf("col", lp.scale_.col).vf("row", lp.scale_.row);
  o.o("scale", s);
  return o;
}

Obj statusObj(const HighsSimplexStatus& st) {
  Obj o;
  o.b("initialised_for_new_lp", st.initialised_for_new_lp)
      .b("is_dualized", st.is_dualized).b("is_permuted", st.is_permuted)
      .b("initialised_for_solve", st.initialised_for_solve)
      .b("has_basis", st.has_basis).b("has_ar_matrix", st.has_ar_matrix)
      .b("has_nla", st.has_nla)
      .b("has_dual_steepest_edge_weights", st.has_dual_steepest_edge_weights)
      .b("has_invert", st.has_invert).b("has_fresh_invert", st.has_fresh_invert)
      .b("has_fresh_rebuild", st.has_fresh_rebuild)
      .b("has_dual_objective_value", st.has_dual_objective_value)
      .b("has_primal_objective_value", st.has_primal_objective_value);
  return o;
}

Obj basisObj(const SimplexBasis& b) {
  Obj o;
  o.vi("basicIndex", b.basicIndex_).vi8("nonbasicFlag", b.nonbasicFlag_)
      .vi8("nonbasicMove", b.nonbasicMove_);
  return o;
}

// Hashes of the solver work arrays: localise a divergence to an iteration.
Obj hashObj(const HEkk& e) {
  const HighsSimplexInfo& in = e.info_;
  Obj o;
  o.h("basicIndex", fnvI(e.basis_.basicIndex_))
      .h("nonbasicMove", fnvI8(e.basis_.nonbasicMove_))
      .h("workCost", fnvF(in.workCost_)).h("workDual", fnvF(in.workDual_))
      .h("workShift", fnvF(in.workShift_)).h("workLower", fnvF(in.workLower_))
      .h("workUpper", fnvF(in.workUpper_)).h("workValue", fnvF(in.workValue_))
      .h("baseValue", fnvF(in.baseValue_)).h("baseLower", fnvF(in.baseLower_))
      .h("baseUpper", fnvF(in.baseUpper_))
      .h("dual_edge_weight", fnvF(e.dual_edge_weight_));
  return o;
}

Obj workArrays(const HEkk& e) {
  const HighsSimplexInfo& in = e.info_;
  Obj o;
  o.vf("workCost", in.workCost_).vf("workDual", in.workDual_)
      .vf("workShift", in.workShift_).vf("workLower", in.workLower_)
      .vf("workUpper", in.workUpper_).vf("workRange", in.workRange_)
      .vf("workValue", in.workValue_).vf("baseLower", in.baseLower_)
      .vf("baseUpper", in.baseUpper_).vf("baseValue", in.baseValue_)
      .vf("dual_edge_weight", e.dual_edge_weight_);
  return o;
}

Obj ekkScalars(const HEkk& e) {
  const HighsSimplexInfo& in = e.info_;
  Obj o;
  o.i("iteration_count", e.iteration_count_).i("update_count", in.update_count)
      .f("updated_dual_objective_value", in.updated_dual_objective_value)
      .f("updated_primal_objective_value", in.updated_primal_objective_value)
      .f("dual_objective_value", in.dual_objective_value)
      .f("primal_objective_value", in.primal_objective_value)
      .i("num_primal_infeasibilities", in.num_primal_infeasibilities)
      .f("sum_primal_infeasibilities", in.sum_primal_infeasibilities)
      .i("num_dual_infeasibilities", in.num_dual_infeasibilities)
      .f("sum_dual_infeasibilities", in.sum_dual_infeasibilities)
      .b("costs_perturbed", in.costs_perturbed).b("costs_shifted", in.costs_shifted)
      .b("bounds_perturbed", in.bounds_perturbed).b("bounds_shifted", in.bounds_shifted)
      .f("cost_perturbation_base", e.cost_perturbation_base_)
      .f("cost_perturbation_max_abs_cost", e.cost_perturbation_max_abs_cost_)
      .f("build_synthetic_tick", e.build_synthetic_tick_)
      .f("total_synthetic_tick", e.total_synthetic_tick_)
      .f("col_aq_density", in.col_aq_density).f("row_ep_density", in.row_ep_density)
      .f("row_ap_density", in.row_ap_density).f("row_DSE_density", in.row_DSE_density)
      .f("col_BFRT_density", in.col_BFRT_density)
      .f("primal_col_density", in.primal_col_density)
      .f("dual_col_density", in.dual_col_density)
      .i("dual_simplex_cleanup_level", e.dual_simplex_cleanup_level_)
      .b("simplex_in_scaled_space", e.simplex_in_scaled_space_);
  return o;
}

void ekkHook(int event, HEkk& e) {
  if (g.snap) {
    if (event == kOracleEkkSolveExit && !g.snap->valid) {
      g.snap->scaled_lp = e.lp_;
      g.snap->basic_index = e.basis_.basicIndex_;
      g.snap->nonbasic_flag = e.basis_.nonbasicFlag_;
      g.snap->iteration_count = e.iteration_count_;
      g.snap->valid = true;
    }
    return;
  }
  if (!g.out) return;
  Out& out = *g.out;
  switch (event) {
    case kOracleEkkSolveEnter: {
      g.solve_call++;
      g.last_iter = e.iteration_count_;
      g.pending_reason = 0;
      Obj o("ekk_solve_enter");
      o.i("solve_call", g.solve_call).i("iteration_count", e.iteration_count_)
          .o("status", statusObj(e.status_)).o("options", optionsObj(*e.options_))
          .f("cost_scale", e.cost_scale_).o("lp", lpObj(e.lp_));
      if (e.status_.has_basis) o.o("basis", basisObj(e.basis_));
      out.put(o);
      break;
    }
    case kOracleEkkSolveInitialised: {
      const HighsSimplexInfo& in = e.info_;
      Obj o("ekk_solve_initialised");
      o.i("solve_call", g.solve_call).o("status", statusObj(e.status_))
          .o("basis", basisObj(e.basis_)).o("scalars", ekkScalars(e))
          .vf("numTotRandomValue", in.numTotRandomValue_)
          .vi("numTotPermutation", in.numTotPermutation_)
          .vi("numColPermutation", in.numColPermutation_)
          .i("simplex_strategy", in.simplex_strategy)
          .i("dual_edge_weight_strategy", in.dual_edge_weight_strategy)
          .i("primal_edge_weight_strategy", in.primal_edge_weight_strategy)
          .i("price_strategy", in.price_strategy)
          .f("dual_simplex_cost_perturbation_multiplier",
             in.dual_simplex_cost_perturbation_multiplier)
          .f("primal_simplex_bound_perturbation_multiplier",
             in.primal_simplex_bound_perturbation_multiplier)
          .f("factor_pivot_threshold", in.factor_pivot_threshold)
          .i("update_limit", in.update_limit)
          .o("arrays", workArrays(e)).mat("ar_matrix", e.ar_matrix_);
      out.put(o);
      break;
    }
    case kOracleEkkComputeFactor: {
      const HFactor& f = HighsOracleAccess::factor(e);
      Obj o("invert");
      o.i("solve_call", g.solve_call).i("iteration_count", e.iteration_count_)
          .i("rank_deficiency", f.rank_deficiency)
          .f("factor_build_synthetic_tick", f.build_synthetic_tick)
          .f("build_synthetic_tick", e.build_synthetic_tick_)
          .f("total_synthetic_tick", e.total_synthetic_tick_)
          .i("basis_matrix_num_el", f.basis_matrix_num_el)
          .i("invert_num_el", f.invert_num_el).i("kernel_dim", f.kernel_dim)
          .i("kernel_num_el", f.kernel_num_el)
          .b("has_invert", e.status_.has_invert)
          .h("h_basicIndex", fnvI(e.basis_.basicIndex_));
      out.put(o);
      break;
    }
    case kOracleEkkSolveExit: {
      Obj o("ekk_solve_exit");
      o.i("solve_call", g.solve_call).i("model_status", (int)e.model_status_)
          .i("exit_algorithm", (int)e.exit_algorithm_)
          .b("solve_bailout", e.solve_bailout_).o("scalars", ekkScalars(e))
          .o("status", statusObj(e.status_)).o("basis", basisObj(e.basis_))
          .o("arrays", workArrays(e));
      out.put(o);
      break;
    }
  }
}

void dualHook(int event, HEkkDual& d) {
  if (!g.out) return;
  Out& out = *g.out;
  const HEkk& e = d.ekk_instance_;
  switch (event) {
    case kOracleSolveStart: {
      g.dual_call++;
      Obj o("dual_start");
      o.i("solve_call", g.solve_call).i("dual_call", g.dual_call)
          .i("solve_phase", d.solve_phase).i("edge_weight_mode", (int)d.edge_weight_mode)
          .b("initial_basis_is_logical", d.initial_basis_is_logical_);
      // interpretDualEdgeWeightStrategy sets this member only on the steepest
      // edge branches; with Devex/Dantzig (e.g. HApp forcing Devex for the
      // unscaled re-solve) it is uninitialised and never read, so omit it.
      if (d.edge_weight_mode == EdgeWeightMode::kSteepestEdge)
        o.b("allow_dual_steepest_edge_to_devex_switch",
            d.allow_dual_steepest_edge_to_devex_switch);
      o.f("Tp", d.Tp).f("Td", d.Td).f("objective_bound", d.objective_bound)
          .b("force_phase2", d.force_phase2).o("scalars", ekkScalars(e))
          .o("arrays", workArrays(e));
      out.put(o);
      break;
    }
    case kOracleRebuildEnter:
      g.pending_reason = d.rebuild_reason;
      break;
    case kOracleRebuild: {
      Obj o("rebuild");
      o.s("algorithm", "dual").i("solve_call", g.solve_call)
          .i("reason", g.pending_reason).i("solve_phase", d.solve_phase)
          .i("rebuild_reason_after", d.rebuild_reason)
          .i("dualInfeasCount", d.dualInfeasCount)
          .i("edge_weight_mode", (int)d.edge_weight_mode)
          .o("scalars", ekkScalars(e)).o("hash", hashObj(e));
      out.put(o);
      g.pending_reason = 0;
      break;
    }
    case kOracleIterate: {
      const bool done = e.iteration_count_ != g.last_iter;
      g.last_iter = e.iteration_count_;
      Obj o(done ? "iter" : "iter_none");
      o.s("algorithm", "dual").i("iteration", e.iteration_count_)
          .i("solve_phase", d.solve_phase).i("rebuild_reason", d.rebuild_reason)
          .i("row_out", d.row_out);
      if (done) {
        std::vector<HighsInt> flips;
        for (HighsInt k = 0; k < d.dualRow.workCount; k++)
          flips.push_back(d.dualRow.workData[k].first);
        o.i("variable_out", d.variable_out).i("variable_in", d.variable_in)
            .i("move_out", d.move_out).f("delta_primal", d.delta_primal)
            .f("theta_dual", d.theta_dual).f("theta_primal", d.theta_primal)
            .f("alpha_col", d.alpha_col).f("alpha_row", d.alpha_row)
            .f("numericalTrouble", d.numericalTrouble)
            .f("computed_edge_weight", d.computed_edge_weight)
            .i("edge_weight_mode", (int)d.edge_weight_mode)
            .i("num_flips", d.dualRow.workCount).vi("flips", flips)
            .i("update_count", e.info_.update_count)
            .f("updated_dual_objective_value", e.info_.updated_dual_objective_value)
            .f("total_synthetic_tick", e.total_synthetic_tick_)
            .f("row_ep_density", e.info_.row_ep_density)
            .f("row_ap_density", e.info_.row_ap_density)
            .f("col_aq_density", e.info_.col_aq_density)
            .f("row_DSE_density", e.info_.row_DSE_density);
        if (d.row_out >= 0 && d.row_out < (HighsInt)e.dual_edge_weight_.size())
          o.f("edge_weight_row_out_after", e.dual_edge_weight_[d.row_out]);
        o.h("h_basicIndex", fnvI(e.basis_.basicIndex_))
            .h("h_workDual", fnvF(e.info_.workDual_))
            .h("h_baseValue", fnvF(e.info_.baseValue_))
            .h("h_dual_edge_weight", fnvF(e.dual_edge_weight_));
      }
      out.put(o);
      break;
    }
    case kOracleCleanup: {
      Obj o("cleanup");
      o.s("algorithm", "dual").i("solve_call", g.solve_call)
          .i("solve_phase", d.solve_phase).o("scalars", ekkScalars(e))
          .o("hash", hashObj(e));
      out.put(o);
      break;
    }
    case kOraclePhaseEnd: {
      Obj o("phase_end");
      o.s("algorithm", "dual").i("solve_call", g.solve_call)
          .i("solve_phase", d.solve_phase).i("rebuild_reason", d.rebuild_reason)
          .i("model_status", (int)e.model_status_)
          .i("iteration", e.iteration_count_);
      out.put(o);
      break;
    }
  }
}

void primalHook(int event, HEkkPrimal& p) {
  if (!g.out) return;
  Out& out = *g.out;
  const HEkk& e = HighsOracleAccess::ekk(p);
  switch (event) {
    case kOracleSolveStart: {
      g.primal_call++;
      g.primal_iterated = false;
      Obj o("primal_start");
      o.i("solve_call", g.solve_call).i("primal_call", g.primal_call)
          .o("state", HighsOracleAccess::primalState(p)).o("scalars", ekkScalars(e))
          .o("arrays", workArrays(e));
      out.put(o);
      break;
    }
    case kOracleRebuildEnter:
      g.pending_reason = HighsOracleAccess::rebuildReason(p);
      break;
    case kOracleRebuild: {
      Obj o("rebuild");
      o.s("algorithm", "primal").i("solve_call", g.solve_call)
          .i("reason", g.pending_reason).o("state", HighsOracleAccess::primalState(p))
          .o("scalars", ekkScalars(e)).o("hash", hashObj(e));
      out.put(o);
      g.pending_reason = 0;
      break;
    }
    case kOracleIterate: {
      const bool done = e.iteration_count_ != g.last_iter;
      g.last_iter = e.iteration_count_;
      if (done) g.primal_iterated = true;
      Obj it = HighsOracleAccess::primalIter(p, g.primal_iterated);
      Obj o(done ? "iter" : "iter_none");
      o.s("algorithm", "primal").i("iteration", e.iteration_count_).o("p", it);
      if (done) {
        o.i("update_count", e.info_.update_count)
            .f("updated_primal_objective_value", e.info_.updated_primal_objective_value)
            .f("total_synthetic_tick", e.total_synthetic_tick_)
            .h("h_basicIndex", fnvI(e.basis_.basicIndex_))
            .h("h_workDual", fnvF(e.info_.workDual_))
            .h("h_baseValue", fnvF(e.info_.baseValue_))
            .h("h_workValue", fnvF(e.info_.workValue_));
      }
      out.put(o);
      break;
    }
    case kOracleCleanup: {
      Obj o("cleanup");
      o.s("algorithm", "primal").i("solve_call", g.solve_call)
          .o("state", HighsOracleAccess::primalState(p)).o("scalars", ekkScalars(e))
          .o("hash", hashObj(e));
      out.put(o);
      break;
    }
    case kOraclePhaseEnd: {
      Obj o("phase_end");
      o.s("algorithm", "primal").i("solve_call", g.solve_call)
          .o("state", HighsOracleAccess::primalState(p))
          .i("model_status", (int)e.model_status_)
          .i("iteration", e.iteration_count_);
      out.put(o);
      break;
    }
  }
}

void installHooks() {
  highsOracleHooks().ekk = ekkHook;
  highsOracleHooks().dual = dualHook;
  highsOracleHooks().primal = primalHook;
}

void removeHooks() {
  highsOracleHooks().ekk = nullptr;
  highsOracleHooks().dual = nullptr;
  highsOracleHooks().primal = nullptr;
}

}  // namespace

void oracleBaseOptions(Highs& highs, int random_seed) {
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("threads", 1);
  highs.setOptionValue("parallel", "off");
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("solver", "simplex");
  highs.setOptionValue("random_seed", random_seed);
}

HighsLp oracleLoadLp(const std::string& path) {
  Highs highs;
  oracleBaseOptions(highs, 0);
  if (highs.readModel(path) == HighsStatus::kError) {
    fprintf(stderr, "highs-oracle: cannot read %s\n", path.c_str());
    exit(2);
  }
  HighsLp lp = highs.getLp();
  lp.integrality_.clear();
  return lp;
}

static bool readRelaxation(Highs& highs, const std::string& path) {
  if (highs.readModel(path) == HighsStatus::kError) return false;
  highs.clearIntegrality();
  return true;
}

EkkSnapshot oracleSnapshot(const std::string& path, HighsInt iteration_limit) {
  EkkSnapshot snap;
  Highs highs;
  oracleBaseOptions(highs, 0);
  highs.setOptionValue("simplex_strategy", 1);
  if (iteration_limit >= 0) highs.setOptionValue("simplex_iteration_limit", iteration_limit);
  if (!readRelaxation(highs, path)) return snap;
  g = TraceState();
  g.snap = &snap;
  installHooks();
  highs.run();
  removeHooks();
  g = TraceState();
  return snap;
}

int runSimplex(const std::string& mps, const std::string& name, int strategy,
               const std::string& path) {
  return runSimplexThreads(mps, name, strategy, 0, path);
}

// threads = 0: the serial oracle (threads=1, parallel=off) and the original
// trace_header. threads > 0 (w5-simplex, w5.cpp): threads=N, parallel=on, and
// the header also records "threads" and "parallel" (FORMAT.md "W5 records").
int runSimplexThreads(const std::string& mps, const std::string& name, int strategy,
                      int threads, const std::string& path) {
  Highs highs;
  const int seed = 0;  // HiGHS default random_seed, set explicitly
  oracleBaseOptions(highs, seed);
  highs.setOptionValue("simplex_strategy", strategy);
  if (threads > 0) {
    highs.setOptionValue("threads", threads);
    highs.setOptionValue("parallel", "on");
  }
  if (!readRelaxation(highs, mps)) {
    fprintf(stderr, "highs-oracle: cannot read %s\n", mps.c_str());
    return 2;
  }
  Out out(path);
  Obj header("trace_header");
  header.s("name", name).i("simplex_strategy", strategy);
  if (threads > 0) header.i("threads", threads).s("parallel", "on");
  header.i("random_seed", seed).s("highs_version", highsVersion())
      .s("highs_githash", highsGithash());
  out.put(header);
  g = TraceState();
  g.out = &out;
  installHooks();
  HighsStatus run_status = highs.run();
  removeHooks();
  g = TraceState();

  const HighsInfo& info = highs.getInfo();
  const HighsSolution& sol = highs.getSolution();
  const HighsBasis& basis = highs.getBasis();
  std::vector<HighsInt> col_status, row_status;
  for (auto s : basis.col_status) col_status.push_back((HighsInt)s);
  for (auto s : basis.row_status) row_status.push_back((HighsInt)s);
  out.put(Obj("result").i("run_status", (int)run_status)
              .i("model_status", (int)highs.getModelStatus())
              .s("model_status_string", highs.modelStatusToString(highs.getModelStatus()))
              .f("objective_function_value", info.objective_function_value)
              .i("simplex_iteration_count", info.simplex_iteration_count)
              .i("primal_solution_status", info.primal_solution_status)
              .i("dual_solution_status", info.dual_solution_status)
              .i("num_primal_infeasibilities", info.num_primal_infeasibilities)
              .i("num_dual_infeasibilities", info.num_dual_infeasibilities)
              .b("basis_valid", basis.valid).vi("col_status", col_status)
              .vi("row_status", row_status).vf("col_value", sol.col_value)
              .vf("col_dual", sol.col_dual).vf("row_value", sol.row_value)
              .vf("row_dual", sol.row_dual));
  return 0;
}
