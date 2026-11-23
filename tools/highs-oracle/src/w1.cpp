// W1 oracle commands: QUASS (w1-qp), getRanging (w1-ranging), getIis (w1-iis)
// and iCrash (w1-icrash), all through the public Highs API on a fresh Highs
// with the W1 base options. Records: FORMAT.md "W1 records".
//
//   highs-oracle w1-qp      <file.mps> <name> <out.jsonl>
//   highs-oracle w1-ranging <file.mps> <name> <out.jsonl>
//   highs-oracle w1-iis     <file.mps> <name> <out.jsonl>
//   highs-oracle w1-icrash  <file.mps> <name> <out.jsonl>
//
// w1-qp needs the third patch, highs-oracle-w1-hooks.patch
// (qpsolver/HighsOracleQpDump.h): it writes the w1_header line, points
// HIGHS_ORACLE_QP_DUMP at the fixture itself (the hook appends), runs, and
// appends w1_result, so the qp_iter / qp_end hook lines land in place.
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "Highs.h"
#include "oracle.h"
#include "qpsolver/HighsOracleQpDump.h"

namespace {

// The W1 base options (tests/port/w1.rs w1_options): output off, serial,
// presolve off, random_seed=0. Unlike oracleBaseOptions, solver is left at
// its default ("choose"), so a QP is routed exactly as HiGHS routes it.
void w1BaseOptions(Highs& highs) {
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("threads", 1);
  highs.setOptionValue("parallel", "off");
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("random_seed", 0);
}

bool w1Overrides(Highs& highs, const std::string& cfg) {
  if (cfg == "default") return true;
  std::stringstream ss(cfg);
  std::string kv;
  while (std::getline(ss, kv, ',')) {
    const size_t eq = kv.find('=');
    if (eq == std::string::npos ||
        highs.setOptionValue(kv.substr(0, eq), kv.substr(eq + 1)) != HighsStatus::kOk) {
      fprintf(stderr, "highs-oracle: bad option %s\n", kv.c_str());
      return false;
    }
  }
  return true;
}

// Fresh Highs with the base options and <cfg>, or exit(64).
void w1Configure(Highs& highs, const std::string& cfg) {
  w1BaseOptions(highs);
  if (!w1Overrides(highs, cfg)) exit(64);
}

// The model as Highs::readModel stores it; integrality dropped.
void w1Read(const std::string& mps, HighsLp& lp, HighsHessian& hessian) {
  Highs reader;
  w1BaseOptions(reader);
  if (reader.readModel(mps) == HighsStatus::kError) {
    fprintf(stderr, "highs-oracle: cannot read %s\n", mps.c_str());
    exit(2);
  }
  lp = reader.getLp();
  lp.integrality_.clear();
  hessian = reader.getModel().hessian_;
}

// HighsLp, encoded as the "lp" object of FORMAT.md (unscaled).
Obj w1LpObj(const HighsLp& lp) {
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

Obj w1Header(const std::string& name, const std::string& overrides, const HighsLp& lp,
             const HighsHessian* hessian) {
  Obj h("w1_header");
  h.s("name", name).s("highs_version", highsVersion()).s("highs_githash", highsGithash())
      .s("overrides", overrides).o("lp", w1LpObj(lp));
  if (hessian) {
    Obj q;
    q.i("dim", hessian->dim_).i("format", (int)hessian->format_)
        .vi("start", hessian->start_).vi("index", hessian->index_)
        .vf("value", hessian->value_);
    h.o("hessian", q);
  }
  return h;
}

Obj w1Result(HighsStatus run_status, const Highs& highs) {
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& sol = highs.getSolution();
  const HighsBasis& basis = highs.getBasis();
  std::vector<HighsInt> col_status, row_status;
  for (auto s : basis.col_status) col_status.push_back((HighsInt)s);
  for (auto s : basis.row_status) row_status.push_back((HighsInt)s);
  Obj r("w1_result");
  r.i("run_status", (int)run_status)
      .i("model_status", (int)highs.getModelStatus())
      .f("objective_function_value", info.objective_function_value)
      .i("simplex_iteration_count", info.simplex_iteration_count)
      .i("qp_iteration_count", info.qp_iteration_count)
      .i("primal_solution_status", info.primal_solution_status)
      .i("dual_solution_status", info.dual_solution_status)
      .i("num_primal_infeasibilities", info.num_primal_infeasibilities)
      .f("max_primal_infeasibility", info.max_primal_infeasibility)
      .f("sum_primal_infeasibilities", info.sum_primal_infeasibilities)
      .i("num_dual_infeasibilities", info.num_dual_infeasibilities)
      .f("max_dual_infeasibility", info.max_dual_infeasibility)
      .f("sum_dual_infeasibilities", info.sum_dual_infeasibilities)
      .b("basis_valid", basis.valid).vi("col_status", col_status)
      .vi("row_status", row_status).vf("col_value", sol.col_value)
      .vf("col_dual", sol.col_dual).vf("row_value", sol.row_value)
      .vf("row_dual", sol.row_dual);
  return r;
}

bool w1Append(const std::string& path, const Obj& o) {
  std::FILE* f = std::fopen(path.c_str(), "ab");
  if (f == nullptr) return false;
  const std::string line = o.str() + "\n";
  const bool ok = std::fwrite(line.data(), 1, line.size(), f) == line.size();
  return std::fclose(f) == 0 && ok;
}

void w1Pass(Highs& highs, const HighsLp& lp, const std::string& name) {
  if (highs.passModel(lp) != HighsStatus::kOk) {
    fprintf(stderr, "highs-oracle: passModel failed for %s\n", name.c_str());
    exit(2);
  }
}

void rangingRec(Out& out, const char* which, const HighsRangingRecord& r) {
  Obj o("ranging");
  o.s("which", which).vf("value", r.value_).vf("objective", r.objective_)
      .vi("in_var", r.in_var_).vi("ou_var", r.ou_var_);
  out.put(o);
}

}  // namespace

// w1-qp: QUASS through Highs::run (passModel + passHessian).
int runW1Qp(const std::string& mps, const std::string& name, const std::string& path) {
  const std::string overrides = "default";
  HighsLp lp;
  HighsHessian hessian;
  w1Read(mps, lp, hessian);
  if (hessian.dim_ == 0) {
    fprintf(stderr, "highs-oracle: %s has no Hessian\n", mps.c_str());
    return 2;
  }
  {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
      fprintf(stderr, "highs-oracle: cannot open %s\n", path.c_str());
      return 2;
    }
    const std::string h = w1Header(name, overrides, lp, &hessian).str() + "\n";
    std::fwrite(h.data(), 1, h.size(), f);
    if (std::fclose(f) != 0) return 2;
  }
  Highs highs;
  w1Configure(highs, overrides);
  w1Pass(highs, lp, name);
  if (highs.passHessian(hessian) != HighsStatus::kOk) {
    fprintf(stderr, "highs-oracle: passHessian failed for %s\n", name.c_str());
    return 2;
  }
  setenv("HIGHS_ORACLE_QP_DUMP", path.c_str(), 1);
  highsOracleQpCall() = 0;
  const HighsStatus run_status = highs.run();
  std::FILE* hook = highsOracleQpFile();  // opens it if QUASS never wrote
  unsetenv("HIGHS_ORACLE_QP_DUMP");
  if (hook == nullptr) {
    fprintf(stderr, "highs-oracle: QP hook file not open (%s)\n", path.c_str());
    return 3;
  }
  // No hook line is legitimate: a QP without constraints is solved by
  // computeStartingPointBounded (a_quass.cpp) and never enters Quass::solve.
  std::fflush(hook);
  if (!w1Append(path, w1Result(run_status, highs))) return 2;
  return 0;
}

// w1-ranging: run (simplex) then getRanging.
int runW1Ranging(const std::string& mps, const std::string& name, const std::string& path) {
  const std::string overrides = "default";
  HighsLp lp;
  HighsHessian hessian;
  w1Read(mps, lp, hessian);
  Out out(path);
  if (!out.ok()) return 2;
  out.put(w1Header(name, overrides, lp, nullptr));
  Highs highs;
  w1Configure(highs, overrides);
  w1Pass(highs, lp, name);
  const HighsStatus run_status = highs.run();
  out.put(w1Result(run_status, highs));
  HighsRanging ranging;
  const HighsStatus status = highs.getRanging(ranging);
  Obj rs("ranging_status");
  rs.i("status", (int)status).b("valid", ranging.valid);
  out.put(rs);
  rangingRec(out, "col_cost_up", ranging.col_cost_up);
  rangingRec(out, "col_cost_dn", ranging.col_cost_dn);
  rangingRec(out, "col_bound_up", ranging.col_bound_up);
  rangingRec(out, "col_bound_dn", ranging.col_bound_dn);
  rangingRec(out, "row_bound_up", ranging.row_bound_up);
  rangingRec(out, "row_bound_dn", ranging.row_bound_dn);
  return 0;
}

// w1-iis: per configuration, a fresh Highs: run then getIis.
int runW1Iis(const std::string& mps, const std::string& name, const std::string& path) {
  static const char* const kCfgs[][2] = {{"light", "iis_strategy=0"},
                                         {"lp", "iis_strategy=2"},
                                         {"ray_irr", "iis_strategy=5"},
                                         {"lp_irr", "iis_strategy=6"}};
  HighsLp lp;
  HighsHessian hessian;
  w1Read(mps, lp, hessian);
  Out out(path);
  if (!out.ok()) return 2;
  out.put(w1Header(name, "default", lp, nullptr));
  for (const auto& c : kCfgs) {
    Highs highs;
    w1Configure(highs, c[1]);
    w1Pass(highs, lp, name);
    const HighsStatus run_status = highs.run();
    const HighsModelStatus model_status = highs.getModelStatus();
    HighsIis iis;
    const HighsStatus iis_status = highs.getIis(iis);
    Obj info;
    info.i("num_lp_solved", iis.info_.num_lp_solved)
        .i("sum_simplex_iteration_counts", iis.info_.sum_simplex_iteration_counts)
        .i("min_simplex_iteration_count", iis.info_.min_simplex_iteration_count)
        .i("max_simplex_iteration_count", iis.info_.max_simplex_iteration_count);
    Obj r("iis_result");
    r.s("cfg", c[0]).s("overrides", c[1]).i("run_status", (int)run_status)
        .i("model_status", (int)model_status).i("get_iis_status", (int)iis_status)
        .b("valid", iis.valid_).i("status", iis.status_).i("strategy", iis.strategy_)
        .vi("col_index", iis.col_index_).vi("row_index", iis.row_index_)
        .vi("col_bound", iis.col_bound_).vi("row_bound", iis.row_bound_)
        .vi("col_status", iis.col_status_).vi("row_status", iis.row_status_)
        .o("info", info);
    out.put(r);
  }
  return 0;
}

// w1-icrash: run with icrash=true (primal simplex crossover).
int runW1Icrash(const std::string& mps, const std::string& name, const std::string& path) {
  const std::string overrides = "icrash=true,simplex_strategy=4";
  HighsLp lp;
  HighsHessian hessian;
  w1Read(mps, lp, hessian);
  Out out(path);
  if (!out.ok()) return 2;
  out.put(w1Header(name, overrides, lp, nullptr));
  Highs highs;
  w1Configure(highs, overrides);
  w1Pass(highs, lp, name);
  const HighsStatus run_status = highs.run();
  const ICrashInfo& ic = highs.getICrashInfo();
  for (const ICrashIterationDetails& d : ic.details) {
    Obj o("icrash_iter");
    o.i("num", d.num).f("weight", d.weight).f("lambda_norm_2", d.lambda_norm_2)
        .f("lp_objective", d.lp_objective).f("quadratic_objective", d.quadratic_objective)
        .f("residual_norm_2", d.residual_norm_2);
    out.put(o);
  }
  Obj ii("icrash_info");
  ii.i("num_iterations", ic.num_iterations).f("final_lp_objective", ic.final_lp_objective)
      .f("final_quadratic_objective", ic.final_quadratic_objective)
      .f("final_residual_norm_2", ic.final_residual_norm_2)
      .f("starting_weight", ic.starting_weight).f("final_weight", ic.final_weight)
      .vf("x_values", ic.x_values);
  out.put(ii);
  out.put(w1Result(run_status, highs));
  return 0;
}
