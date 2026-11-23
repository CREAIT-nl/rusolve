// G1 oracle commands: the Highs API around run (user scaling, infinite costs,
// semi-variables, setSolution, multiple objectives, presolve / postsolve /
// crossover, standard form, rays, ill-conditioning, model modification and
// the improving-solution file), all through the public Highs API on a fresh
// Highs with the G1 base options. Records: FORMAT.md "G1 records".
//
//   highs-oracle g1-<x> <model> <name> <scenario> <out.jsonl>
//
// with <x> one of user-scale, inf-cost, semi, set-solution,
// multi-objective, presolve-api, standard-form, rays, ill-conditioning,
// modify, mipfile, clique-normalise (clqn.cpp, 1.16 only). No hook patch is
// needed. Each command writes the
// g1_header, then one g1_call per API call and one g1_result after every
// run / postsolve / crossover, in the order tests/port/g1.rs replays them.
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Highs.h"
#include "oracle.h"

namespace {

// The G1 base options (tests/port/g1.rs g1_options): output off, serial,
// random_seed=0. Presolve is left at its default unless overridden.
void g1BaseOptions(Highs& highs) {
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("threads", 1);
  highs.setOptionValue("parallel", "off");
  highs.setOptionValue("random_seed", 0);
}

// Applies "k=v,k=v" (or "default"). The mip_improving_solution_file value
// is a bare file name; it is written under file_dir.
bool g1Overrides(Highs& highs, const std::string& cfg, const std::string& file_dir) {
  if (cfg == "default") return true;
  std::stringstream ss(cfg);
  std::string kv;
  while (std::getline(ss, kv, ',')) {
    const size_t eq = kv.find('=');
    if (eq == std::string::npos) {
      fprintf(stderr, "highs-oracle: bad option %s\n", kv.c_str());
      return false;
    }
    const std::string k = kv.substr(0, eq);
    std::string v = kv.substr(eq + 1);
    if (k == "mip_improving_solution_file") v = file_dir + "/" + v;
    if (highs.setOptionValue(k, v) != HighsStatus::kOk) {
      fprintf(stderr, "highs-oracle: bad option %s\n", kv.c_str());
      return false;
    }
  }
  return true;
}

// The model as Highs::readModel stores it, reduced to the fields of the
// FORMAT.md "lp" object plus the integrality (no names: the Rust replay
// passes exactly these fields to passModel).
HighsLp g1Read(const std::string& path) {
  Highs reader;
  g1BaseOptions(reader);
  if (reader.readModel(path) == HighsStatus::kError) {
    fprintf(stderr, "highs-oracle: cannot read %s\n", path.c_str());
    exit(2);
  }
  const HighsLp& r = reader.getLp();
  HighsLp lp;
  lp.num_col_ = r.num_col_;
  lp.num_row_ = r.num_row_;
  lp.col_cost_ = r.col_cost_;
  lp.col_lower_ = r.col_lower_;
  lp.col_upper_ = r.col_upper_;
  lp.row_lower_ = r.row_lower_;
  lp.row_upper_ = r.row_upper_;
  lp.a_matrix_ = r.a_matrix_;
  lp.a_matrix_.ensureColwise();
  lp.sense_ = r.sense_;
  lp.offset_ = r.offset_;
  lp.scale_ = r.scale_;
  lp.is_scaled_ = r.is_scaled_;
  lp.is_moved_ = r.is_moved_;
  lp.integrality_ = r.integrality_;
  return lp;
}

// HighsLp, encoded as the "lp" object of FORMAT.md (column-wise copy).
Obj g1LpObj(const HighsLp& in) {
  HighsSparseMatrix a = in.a_matrix_;
  a.ensureColwise();
  Obj o;
  o.i("num_col", in.num_col_).i("num_row", in.num_row_).i("sense", (int)in.sense_)
      .f("offset", in.offset_).vf("col_cost", in.col_cost_)
      .vf("col_lower", in.col_lower_).vf("col_upper", in.col_upper_)
      .vf("row_lower", in.row_lower_).vf("row_upper", in.row_upper_)
      .mat("a", a).b("is_scaled", in.is_scaled_).b("is_moved", in.is_moved_);
  Obj s;
  s.i("strategy", in.scale_.strategy).b("has_scaling", in.scale_.has_scaling)
      .i("num_col", in.scale_.num_col).i("num_row", in.scale_.num_row)
      .f("cost", in.scale_.cost).vf("col", in.scale_.col).vf("row", in.scale_.row);
  o.o("scale", s);
  return o;
}

// "lp" plus "integrality" (only when non-empty) on o.
void g1AddLp(Obj& o, const HighsLp& lp) {
  o.o("lp", g1LpObj(lp));
  if (!lp.integrality_.empty()) {
    std::vector<HighsInt> t;
    for (HighsVarType v : lp.integrality_) t.push_back((HighsInt)v);
    o.vi("integrality", t);
  }
}

Obj g1LpWithIntegrality(const HighsLp& lp) {
  Obj o;
  g1AddLp(o, lp);
  return o;
}

std::vector<HighsInt> statusInts(const std::vector<HighsBasisStatus>& v) {
  std::vector<HighsInt> r;
  for (HighsBasisStatus s : v) r.push_back((HighsInt)s);
  return r;
}

// A HighsSolution exactly as tests/port/g1.rs solution_from builds it from
// the record fields that are present.
HighsSolution g1Solution(const std::vector<double>* col_value,
                         const std::vector<double>* row_value,
                         const std::vector<double>* col_dual,
                         const std::vector<double>* row_dual) {
  HighsSolution s;
  if (col_value) {
    s.col_value = *col_value;
    s.value_valid = true;
  }
  if (row_value) s.row_value = *row_value;
  if (col_dual) {
    s.col_dual = *col_dual;
    s.row_dual = *row_dual;
    s.dual_valid = true;
  }
  return s;
}

// The solution fields of a g1_call (present iff the vectors were set).
void g1AddSolution(Obj& o, const HighsSolution& s, bool row_value, bool dual) {
  o.vf("col_value", s.col_value);
  if (row_value) o.vf("row_value", s.row_value);
  if (dual) o.vf("col_dual", s.col_dual).vf("row_dual", s.row_dual);
}

class G1 {
 public:
  G1(const std::string& model, const std::string& name, const std::string& overrides,
     const std::string& path)
      : name_(name), overrides_(overrides), out_(path) {
    lp0_ = g1Read(model);
    char tmpl[] = "/tmp/highs-oracle-g1-XXXXXX";
    if (overrides.find("mip_improving_solution_file") != std::string::npos) {
      if (mkdtemp(tmpl) == nullptr) {
        fprintf(stderr, "highs-oracle: mkdtemp failed\n");
        exit(2);
      }
      tmp_dir_ = tmpl;
    }
    configure(highs_, overrides);
    Obj h("g1_header");
    h.s("name", name).s("highs_version", highsVersion()).s("highs_githash", highsGithash())
        .s("overrides", overrides);
    g1AddLp(h, lp0_);
    out_.put(h);
    if (highs_.passModel(lp0_) != HighsStatus::kOk) {
      fprintf(stderr, "highs-oracle: passModel failed for %s\n", name.c_str());
      exit(2);
    }
  }
  ~G1() {
    if (!tmp_dir_.empty()) {
      std::remove(improvingFile().c_str());
      rmdir(tmp_dir_.c_str());
    }
  }

  void configure(Highs& h, const std::string& overrides) {
    g1BaseOptions(h);
    if (!g1Overrides(h, overrides, tmp_dir_)) exit(64);
  }

  // A fresh g1_call record (call, op, status); the caller adds the op fields.
  Obj call(const char* op, HighsStatus status) {
    Obj c("g1_call");
    c.i("call", call_++).s("op", op).i("status", (int)status);
    return c;
  }
  void put(const Obj& o) { out_.put(o); }

  // g1_result after run / postsolve / crossover: the incumbent state.
  void result() {
    const HighsInfo& info = highs_.getInfo();
    const HighsSolution& sol = highs_.getSolution();
    const HighsBasis& basis = highs_.getBasis();
    const HighsLp& lp = highs_.getLp();
    Obj r("g1_result");
    r.i("model_status", (int)highs_.getModelStatus())
        .f("objective_function_value", info.objective_function_value)
        .i("primal_solution_status", info.primal_solution_status)
        .i("dual_solution_status", info.dual_solution_status)
        .i("simplex_iteration_count", info.simplex_iteration_count)
        .f("max_primal_infeasibility", info.max_primal_infeasibility)
        .f("max_dual_infeasibility", info.max_dual_infeasibility)
        .b("value_valid", sol.value_valid)
        .b("dual_valid", sol.dual_valid)
        .b("basis_valid", basis.valid)
        .vf("col_value", sol.col_value).vf("row_value", sol.row_value)
        .vf("col_dual", sol.col_dual).vf("row_dual", sol.row_dual);
    if (lp.isMip()) {
      r.i("mip_node_count", info.mip_node_count)
          .f("mip_dual_bound", info.mip_dual_bound)
          .f("mip_gap", info.mip_gap)
          .f("max_integrality_violation", info.max_integrality_violation);
    }
    if (basis.valid) {
      r.vi("col_status", statusInts(basis.col_status))
          .vi("row_status", statusInts(basis.row_status));
    }
    g1AddLp(r, lp);
    out_.put(r);
  }

  void run() {
    put(call("run", highs_.run()));
    result();
  }

  std::string improvingFile() const { return tmp_dir_ + "/g1_improve.sol"; }

  Highs highs_;
  HighsLp lp0_;
  std::string name_, overrides_, tmp_dir_;

 private:
  Out out_;
  long long call_ = 0;
};

// The optimal solution of the header model, solved on a separate Highs with
// the base options (and `extra` overrides).
HighsSolution g1AuxSolve(const HighsLp& lp, const std::string& extra) {
  Highs aux;
  g1BaseOptions(aux);
  if (!g1Overrides(aux, extra, "")) exit(64);
  if (aux.passModel(lp) == HighsStatus::kError) exit(2);
  aux.run();
  return aux.getSolution();
}

// Scenario -> overrides for the commands whose scenarios only set options.
struct Scn {
  const char* name;
  const char* overrides;
};

bool lookup(const Scn* table, const std::string& scenario, std::string& overrides) {
  for (const Scn* s = table; s->name; s++) {
    if (scenario == s->name) {
      overrides = s->overrides;
      return true;
    }
  }
  fprintf(stderr, "highs-oracle: unknown scenario %s\n", scenario.c_str());
  return false;
}

// ------------------------------------------------------------ commands

// (a), (b), (c): options only, then run.
int runOnly(const Scn* table, const std::string& model, const std::string& name,
            const std::string& scenario, const std::string& path) {
  std::string ov;
  if (!lookup(table, scenario, ov)) return 64;
  G1 g(model, name, ov, path);
  g.run();
  return 0;
}

const Scn kUserScale[] = {{"o3", "user_objective_scale=3"},
                          {"b-2", "user_bound_scale=-2"},
                          {"o-4b5", "user_objective_scale=-4,user_bound_scale=5"},
                          {"o40", "user_objective_scale=40"},
                          {"o3b2", "user_objective_scale=3,user_bound_scale=2"},
                          {nullptr, nullptr}};
const Scn kInfCost[] = {{"default", "default"}, {"ic1e15", "infinite_cost=1e15"}, {nullptr, nullptr}};
const Scn kSemi[] = {{"default", "default"}, {"relax", "solve_relaxation=true"}, {nullptr, nullptr}};

// (d) setSolution: the values come from a first solve on a separate Highs.
int setSolutionCmd(const std::string& model, const std::string& name,
                   const std::string& scenario, const std::string& path) {
  const std::string ov = scenario == "cold" ? "use_warm_start=false" : "default";
  if (scenario != "full" && scenario != "full-dual" && scenario != "sparse" &&
      scenario != "cold" && scenario != "mip-full" && scenario != "mip-partial" &&
      scenario != "mip-fractional") {
    fprintf(stderr, "highs-oracle: unknown scenario %s\n", scenario.c_str());
    return 64;
  }
  G1 g(model, name, ov, path);
  const HighsLp& lp = g.lp0_;
  const HighsSolution first = g1AuxSolve(lp, "default");
  if (!first.value_valid) {
    fprintf(stderr, "highs-oracle: no first solution for %s\n", name.c_str());
    return 3;
  }
  if (scenario == "full" || scenario == "cold" || scenario == "mip-full" ||
      scenario == "mip-fractional") {
    std::vector<double> x = first.col_value;
    if (scenario == "mip-fractional") {
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        if (lp.integrality_[iCol] == HighsVarType::kInteger &&
            x[iCol] + 0.5 <= lp.col_upper_[iCol]) {
          x[iCol] += 0.5;
          break;
        }
      }
    }
    const HighsSolution s = g1Solution(&x, nullptr, nullptr, nullptr);
    Obj c = g.call("set_solution", g.highs_.setSolution(s));
    g1AddSolution(c, s, false, false);
    g.put(c);
  } else if (scenario == "full-dual") {
    const HighsSolution s =
        g1Solution(&first.col_value, &first.row_value, &first.col_dual, &first.row_dual);
    Obj c = g.call("set_solution", g.highs_.setSolution(s));
    g1AddSolution(c, s, true, true);
    g.put(c);
  } else {
    // sparse: every other column; mip-partial: every other integer column
    std::vector<HighsInt> index;
    std::vector<double> value;
    HighsInt k = 0;
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      if (scenario == "mip-partial" && lp.integrality_[iCol] != HighsVarType::kInteger)
        continue;
      if (k++ % 2) continue;
      index.push_back(iCol);
      value.push_back(first.col_value[iCol]);
    }
    Obj c = g.call("set_solution_sparse",
                   g.highs_.setSolution((HighsInt)index.size(), index.data(), value.data()));
    c.vi("index", index).vf("value", value);
    g.put(c);
  }
  g.run();
  return 0;
}

// (e) the objectives of each multi-objective instance.
std::vector<HighsLinearObjective> g1Objectives(const std::string& name, const HighsLp& lp) {
  auto mk = [](double weight, double offset, std::vector<double> coefficients,
               double abs_tolerance, double rel_tolerance, HighsInt priority) {
    HighsLinearObjective o;
    o.weight = weight;
    o.offset = offset;
    o.coefficients = coefficients;
    o.abs_tolerance = abs_tolerance;
    o.rel_tolerance = rel_tolerance;
    o.priority = priority;
    return o;
  };
  std::vector<HighsLinearObjective> v;
  if (name == "multiobj_lp") {
    v.push_back(mk(1, 0, {1, 2, 0, 1}, 1e-6, 0, 2));
    v.push_back(mk(-1, 1, {0, 1, 1, 0}, 0, 1e-4, 1));
    v.push_back(mk(0.5, -2, {1, 0, 0, -1}, 1e-6, 1e-6, 0));
  } else if (name == "multiobj_mip") {
    v.push_back(mk(1, 0, {-3, -2, -1}, 1e-6, 0, 1));
    v.push_back(mk(1, 0.5, {0, 1, -1}, 0, 1e-4, 0));
  } else {
    // afiro: the model cost, then all ones
    v.push_back(mk(1, 0, lp.col_cost_, 1e-6, 0, 1));
    v.push_back(mk(-0.5, 0, std::vector<double>(lp.num_col_, 1.0), 0, 1e-4, 0));
  }
  return v;
}

int multiObjectiveCmd(const std::string& model, const std::string& name,
                      const std::string& scenario, const std::string& path) {
  static const Scn table[] = {{"blend", "blend_multi_objectives=true"},
                              {"lex", "blend_multi_objectives=false"},
                              {nullptr, nullptr}};
  std::string ov;
  if (!lookup(table, scenario, ov)) return 64;
  G1 g(model, name, ov, path);
  const std::vector<HighsLinearObjective> objs = g1Objectives(name, g.lp0_);
  {
    Obj c = g.call("pass_linear_objectives",
                   g.highs_.passLinearObjectives((HighsInt)objs.size(), objs.data()));
    std::string list = "[";
    for (size_t i = 0; i < objs.size(); i++) {
      Obj o;
      o.f("weight", objs[i].weight).f("offset", objs[i].offset)
          .vf("coefficients", objs[i].coefficients)
          .f("abs_tolerance", objs[i].abs_tolerance)
          .f("rel_tolerance", objs[i].rel_tolerance).i("priority", objs[i].priority);
      if (i) list += ",";
      list += o.str();
    }
    list += "]";
    c.raw("objectives", list);
    g.put(c);
  }
  g.run();
  g.put(g.call("clear_linear_objectives", g.highs_.clearLinearObjectives()));
  g.run();
  return 0;
}

// The presolve call: model_presolve_status and, when presolve produced one,
// the presolved model.
void presolveCall(G1& g) {
  const HighsStatus s = g.highs_.presolve();
  const HighsPresolveStatus ps = g.highs_.getModelPresolveStatus();
  Obj c = g.call("presolve", s);
  c.i("model_presolve_status", (int)ps);
  if (s != HighsStatus::kError &&
      (ps == HighsPresolveStatus::kNotReduced || ps == HighsPresolveStatus::kReduced ||
       ps == HighsPresolveStatus::kReducedToEmpty || ps == HighsPresolveStatus::kTimeout))
    c.o("presolved", g1LpWithIntegrality(g.highs_.getPresolvedLp()));
  g.put(c);
}

// (f) presolve / postsolve / crossover.
int presolveApiCmd(const std::string& model, const std::string& name,
                   const std::string& scenario, const std::string& path) {
  if (scenario == "presolve") {
    G1 g(model, name, "default", path);
    presolveCall(g);
    g.run();
  } else if (scenario == "postsolve") {
    G1 g(model, name, "default", path);
    presolveCall(g);
    // Solve the presolved model on a second Highs, presolve off
    Highs aux;
    g1BaseOptions(aux);
    aux.setOptionValue("presolve", "off");
    if (aux.passModel(g.highs_.getPresolvedLp()) == HighsStatus::kError) return 2;
    aux.run();
    const HighsSolution& a = aux.getSolution();
    const bool mip = g.lp0_.isMip();
    const bool dual = !mip && a.dual_valid;
    const HighsSolution s = g1Solution(&a.col_value, &a.row_value, dual ? &a.col_dual : nullptr,
                                       dual ? &a.row_dual : nullptr);
    Obj c = g.call("postsolve", g.highs_.postsolve(s));
    g1AddSolution(c, s, true, dual);
    g.put(c);
    g.result();
  } else if (scenario == "crossover" || scenario == "crossover-dual") {
    // crossover passes the IPM primal values only; crossover-dual adds the
    // duals (IPX CrossoverFromStartingPoint then rejects the HiGHS-signed
    // duals: kError)
    G1 g(model, name, "solver=ipm,run_crossover=off", path);
    g.run();
    const HighsSolution& r = g.highs_.getSolution();
    const bool dual = scenario == "crossover-dual";
    const HighsSolution s = g1Solution(&r.col_value, &r.row_value, dual ? &r.col_dual : nullptr,
                                       dual ? &r.row_dual : nullptr);
    Obj c = g.call("crossover", g.highs_.crossover(s));
    g1AddSolution(c, s, true, dual);
    g.put(c);
    g.result();
  } else {
    fprintf(stderr, "highs-oracle: unknown scenario %s\n", scenario.c_str());
    return 64;
  }
  return 0;
}

// (g) getStandardFormLp: sizes first (null arrays), then the arrays.
void standardFormCall(G1& g) {
  HighsInt nc = 0, nr = 0, nz = 0;
  double off = 0;
  HighsStatus s = g.highs_.getStandardFormLp(nc, nr, nz, off);
  std::vector<double> cost, rhs, value;
  std::vector<HighsInt> start, index;
  if (s != HighsStatus::kError) {
    cost.assign(nc, 0);
    rhs.assign(nr, 0);
    start.assign(nc + 1, 0);
    index.assign(nz, 0);
    value.assign(nz, 0);
    s = g.highs_.getStandardFormLp(nc, nr, nz, off, cost.data(), rhs.data(), start.data(),
                                   index.data(), value.data());
  }
  Obj c = g.call("get_standard_form_lp", s);
  c.i("num_col", nc).i("num_row", nr).i("num_nz", nz).f("offset", off).vf("cost", cost)
      .vf("rhs", rhs).vi("start", start).vi("index", index).vf("value", value);
  g.put(c);
}

int standardFormCmd(const std::string& model, const std::string& name,
                    const std::string& scenario, const std::string& path) {
  if (scenario != "default") {
    fprintf(stderr, "highs-oracle: unknown scenario %s\n", scenario.c_str());
    return 64;
  }
  G1 g(model, name, "default", path);
  standardFormCall(g);
  g.run();
  standardFormCall(g);
  return 0;
}

// (h) rays and the fixed LP.
void rayCall(G1& g, const char* op, HighsInt dim) {
  bool has = false;
  std::vector<double> v(dim, 0.0);
  HighsStatus s;
  if (std::string(op) == "get_dual_ray")
    s = g.highs_.getDualRay(has, v.data());
  else if (std::string(op) == "get_primal_ray")
    s = g.highs_.getPrimalRay(has, v.data());
  else
    s = g.highs_.getDualUnboundednessDirection(has, v.data());
  Obj c = g.call(op, s);
  c.i("dim", dim).b("has", has).vf("value", v);
  g.put(c);
}

int raysCmd(const std::string& model, const std::string& name, const std::string& scenario,
            const std::string& path) {
  if (scenario == "simplex" || scenario == "presolve") {
    G1 g(model, name, scenario == "simplex" ? "presolve=off" : "presolve=on", path);
    g.run();
    rayCall(g, "get_dual_ray", g.lp0_.num_row_);
    rayCall(g, "get_dual_unboundedness_direction", g.lp0_.num_col_);
    rayCall(g, "get_primal_ray", g.lp0_.num_col_);
  } else if (scenario == "fixed") {
    G1 g(model, name, "default", path);
    g.run();
    HighsLp fixed;
    const HighsStatus s = g.highs_.getFixedLp(fixed);
    Obj c = g.call("get_fixed_lp", s);
    if (s != HighsStatus::kError) c.o("fixed", g1LpWithIntegrality(fixed));
    g.put(c);
  } else {
    fprintf(stderr, "highs-oracle: unknown scenario %s\n", scenario.c_str());
    return 64;
  }
  return 0;
}

// (i) getIllConditioning after run.
int illConditioningCmd(const std::string& model, const std::string& name,
                       const std::string& scenario, const std::string& path) {
  if (scenario != "m0c" && scenario != "m0v" && scenario != "m1c" && scenario != "m1v") {
    fprintf(stderr, "highs-oracle: unknown scenario %s\n", scenario.c_str());
    return 64;
  }
  const HighsInt method = scenario[1] == '1' ? 1 : 0;
  const bool constraint = scenario[2] == 'c';
  const double bound = 1e-4;
  G1 g(model, name, "default", path);
  g.run();
  HighsIllConditioning ic;
  const HighsStatus s = g.highs_.getIllConditioning(ic, constraint, method, bound);
  std::vector<HighsInt> index;
  std::vector<double> multiplier;
  for (const HighsIllConditioningRecord& r : ic.record) {
    index.push_back(r.index);
    multiplier.push_back(r.multiplier);
  }
  Obj c = g.call("get_ill_conditioning", s);
  c.b("constraint", constraint).i("method", method).f("ill_conditioning_bound", bound)
      .vi("index", index).vf("multiplier", multiplier);
  g.put(c);
  return 0;
}

// (j) model modification, each change followed by a (warm-started) run.
int modifyCmd(const std::string& model, const std::string& name, const std::string& scenario,
              const std::string& path) {
  G1 g(model, name, "default", path);
  const HighsLp& lp = g.lp0_;
  const std::vector<HighsInt>& start = lp.a_matrix_.start_;
  const std::vector<HighsInt>& aindex = lp.a_matrix_.index_;
  const std::vector<double>& avalue = lp.a_matrix_.value_;
  auto coeff = [&](HighsInt row, HighsInt col, double value) {
    Obj c = g.call("change_coeff", g.highs_.changeCoeff(row, col, value));
    c.i("row", row).i("col", col).f("value", value);
    g.put(c);
    g.run();
  };
  if (scenario == "coeff") {
    // An existing nonzero (doubled), a new nonzero, an existing one to zero
    coeff(aindex[start[0]], 0, 2 * avalue[start[0]]);
    HighsInt new_row = 0;
    for (HighsInt iEl = start[1]; iEl < start[2]; iEl++) {
      if (aindex[iEl] != new_row) break;
      new_row++;
    }
    coeff(new_row, 1, 0.5);
    coeff(aindex[start[2]], 2, 0.0);
  } else if (scenario == "scale") {
    auto scale = [&](const char* op, HighsInt ix, double value) {
      const bool col = std::string(op) == "scale_col";
      Obj c = g.call(op, col ? g.highs_.scaleCol(ix, value) : g.highs_.scaleRow(ix, value));
      c.i(col ? "col" : "row", ix).f("value", value);
      g.put(c);
      g.run();
    };
    scale("scale_col", 0, 2.0);
    scale("scale_col", 1, -0.5);
    scale("scale_row", 0, -2.0);
  } else if (scenario == "add") {
    {
      const std::vector<HighsInt> index = {0, 1, 2};
      const std::vector<double> value = {1.0, -1.0, 2.0};
      const double cost = 1.0, lower = 0.0, upper = kHighsInf;
      Obj c = g.call("add_col", g.highs_.addCol(cost, lower, upper, 3, index.data(),
                                                value.data()));
      c.f("cost", cost).f("lower", lower).f("upper", upper).vi("index", index)
          .vf("value", value);
      g.put(c);
      g.run();
    }
    {
      const std::vector<HighsInt> index = {0, 1, 2};
      const std::vector<double> value = {1.0, 1.0, 1.0};
      const double lower = -kHighsInf, upper = 1000.0;
      Obj c = g.call("add_row", g.highs_.addRow(lower, upper, 3, index.data(), value.data()));
      c.f("lower", lower).f("upper", upper).vi("index", index).vf("value", value);
      g.put(c);
      g.run();
    }
  } else {
    fprintf(stderr, "highs-oracle: unknown scenario %s\n", scenario.c_str());
    return 64;
  }
  return 0;
}

// (k) mip_improving_solution_file: run, the file text, the saved solutions.
int mipfileCmd(const std::string& model, const std::string& name, const std::string& scenario,
               const std::string& path) {
  static const Scn table[] = {
      {"dense", "mip_improving_solution_file=g1_improve.sol,mip_improving_solution_save=true"},
      {"sparse",
       "mip_improving_solution_file=g1_improve.sol,mip_improving_solution_save=true,"
       "mip_improving_solution_report_sparse=true"},
      {nullptr, nullptr}};
  std::string ov;
  if (!lookup(table, scenario, ov)) return 64;
  G1 g(model, name, ov, path);
  g.run();
  {
    std::ifstream in(g.improvingFile(), std::ios::binary);
    std::stringstream text;
    text << in.rdbuf();
    Obj f("g1_file");
    f.s("text", text.str());
    g.put(f);
  }
  std::vector<double> objectives;
  for (const HighsObjectiveSolution& s : g.highs_.getSavedMipSolutions())
    objectives.push_back(s.objective);
  Obj c = g.call("get_saved_mip_solutions", HighsStatus::kOk);
  c.vf("objectives", objectives);
  g.put(c);
  return 0;
}

}  // namespace

int runG1(const std::string& cmd, const std::string& model, const std::string& name,
          const std::string& scenario, const std::string& out) {
  if (cmd == "g1-user-scale") return runOnly(kUserScale, model, name, scenario, out);
  if (cmd == "g1-inf-cost") return runOnly(kInfCost, model, name, scenario, out);
  if (cmd == "g1-semi") return runOnly(kSemi, model, name, scenario, out);
  if (cmd == "g1-set-solution") return setSolutionCmd(model, name, scenario, out);
  if (cmd == "g1-multi-objective") return multiObjectiveCmd(model, name, scenario, out);
  if (cmd == "g1-presolve-api") return presolveApiCmd(model, name, scenario, out);
  if (cmd == "g1-standard-form") return standardFormCmd(model, name, scenario, out);
  if (cmd == "g1-rays") return raysCmd(model, name, scenario, out);
  if (cmd == "g1-ill-conditioning") return illConditioningCmd(model, name, scenario, out);
  if (cmd == "g1-modify") return modifyCmd(model, name, scenario, out);
  if (cmd == "g1-mipfile") return mipfileCmd(model, name, scenario, out);
  if (cmd == "g1-clique-normalise") return runClqn(model, name, scenario, out);
  return -1;
}

HighsLp oracleG1Read(const std::string& path) { return g1Read(path); }

void oracleG1AddLp(Obj& o, const HighsLp& lp) { g1AddLp(o, lp); }
