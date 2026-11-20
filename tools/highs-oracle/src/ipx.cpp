// IPX fixtures (ipx-model, plus the linear-algebra component fixtures): IPX of HiGHS 1.15.1 (highs/ipm/ipx, (C)
// 2018 ERGO-Code, distributed with HiGHS under the MIT license) driven
// through its classes, below the LpSolver.
//
//   ipx-model  <mps> <name> <cfg> <out>
//       params, input (the Model::Load arguments exactly as IpxWrapper's
//       fillInIpxData builds them), model (every Model member after Load:
//       computational form AI/AIt/b/c/lb/ub, scaling factors, flipped and
//       boxed variables, user model after scaling, GetInfo) and transforms
//       (PresolveStartingPoint, PresolveIPMStartingPoint,
//       PostsolveInteriorSolution, EvaluateInteriorSolution of fixed
//       pseudo-random vectors).
//   ipx-linalg <mps> <name> <cfg> <out>
//       params, input, then the linear-algebra components on iterates produced
//       by the real C++ IPM (the IPM itself is not replayed):
//       iterate (state after IPM::StartingPoint / the initial IPM / a few
//       main IPM iterations), kkt_diag (KKTSolverDiag Factorize + Solve),
//       starting_basis (StartingBasis: crash basis, maxvolume-free repair,
//       LU factors), basis_ops (SolveDense, MinSingularValue,
//       SymbolicInvert, RandomPermute), kkt_basis (KKTSolverBasis
//       Factorize + Solve: DropPrimal/DropDual, maxvolume, basis updates),
//       basis_load (Basis::Load of a later basis), postsolve (Postprocess,
//       DropToComplementarity, Model postsolve/evaluate of the iterate and
//       the basic solution).
//
// <cfg> is a comma separated list of ipx::Parameters overrides
// (key=value, or "default"). Every record carries its inputs (or the
// state the Rust replay must start from), so the replay needs no MPS reader
// and no IPM.
//
// Private members are read with `#define private public` (this translation
// unit only, after every standard and HiGHS header the IPX headers include).
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <ostream>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <valarray>
#include <vector>

#include <unistd.h>

#include "HConfig.h"
#include "HighsExternalApi.h"
#include "ipm/basiclu/basiclu.h"
#include "io/HighsIO.h"
#include "lp_data/HighsCallback.h"
#include "parallel/HighsParallel.h"
#include "oracle.h"
#include "util/HighsInt.h"
#include "util/HighsRandom.h"

#define private public
#define protected public
#include "ipm/ipx/basiclu_kernel.h"
#include "ipm/ipx/basiclu_wrapper.h"
#include "ipm/ipx/basis.h"
#include "ipm/ipx/control.h"
#include "ipm/ipx/forrest_tomlin.h"
#include "ipm/ipx/ipm.h"
#include "ipm/ipx/iterate.h"
#include "ipm/ipx/kkt_solver_basis.h"
#include "ipm/ipx/kkt_solver_diag.h"
#include "ipm/ipx/model.h"
#include "ipm/ipx/starting_basis.h"
#include "ipm/ipx/utils.h"
#undef private
#undef protected
// Only present in the oracle copy (highs-oracle-ipx-hooks.patch).
#include "ipm/ipx/HighsOracleIpxDump.h"

// IpxWrapper.cpp (not included: IpxWrapper.h pulls in the IPX headers
// before the access override above).
void fillInIpxData(const HighsLp& lp, ipx::Int& num_col, ipx::Int& num_row,
                   double& offset, std::vector<double>& obj,
                   std::vector<double>& col_lb, std::vector<double>& col_ub,
                   std::vector<ipx::Int>& Ap, std::vector<ipx::Int>& Ai,
                   std::vector<double>& Ax, std::vector<double>& rhs,
                   std::vector<char>& constraint_type);

namespace {

using ipx::Int;
using ipx::Vector;

// Basis::Load of HiGHS 1.15.1 (ipm/ipx/basis.cc), which 1.16 removed (the
// 1.16 library no longer calls it). The basis_load record keeps the same
// operation on both oracle versions, so this is a verbatim copy of the 1.15.1
// body over the members the access override above exposes.
Int loadBasis(ipx::Basis& b, const int* basic_status) {
#if HIGHS_VERSION_MAJOR == 1 && HIGHS_VERSION_MINOR < 16
  return b.Load(basic_status);
#else
  const Int m = b.model_.rows();
  const Int n = b.model_.cols();

  // Change member variables only when basis is valid.
  std::vector<Int> basis, map2basis(n + m);
  Int p = 0;
  for (Int j = 0; j < n + m; j++) {
    switch (basic_status[j]) {
      case ipx::Basis::NONBASIC_FIXED:
        map2basis[j] = -2;
        break;
      case ipx::Basis::NONBASIC:
        map2basis[j] = -1;
        break;
      case ipx::Basis::BASIC:
        basis.push_back(j);
        map2basis[j] = p++;
        break;
      case ipx::Basis::BASIC_FREE:
        basis.push_back(j);
        map2basis[j] = p++ + m;
        break;
      default:
        return IPX_ERROR_invalid_basis;
    }
  }
  if (p != m) return IPX_ERROR_invalid_basis;

  std::copy(basis.begin(), basis.end(), b.basis_.begin());
  std::copy(map2basis.begin(), map2basis.end(), b.map2basis_.begin());
  return b.Factorize();
#endif
}

Obj& vv(Obj& o, const char* k, const Vector& v) {
  std::vector<double> c(std::begin(v), std::end(v));
  return o.vf(k, c);
}

Obj matObj(const ipx::SparseMatrix& a) {
  Obj o;
  const Int nz = a.entries();
  std::vector<Int> colptr(a.colptr(), a.colptr() + a.cols() + 1);
  std::vector<Int> rowidx(a.rowidx(), a.rowidx() + nz);
  std::vector<double> values(a.values(), a.values() + nz);
  o.i("rows", a.rows()).i("cols", a.cols()).vi("colptr", colptr)
      .vi("rowidx", rowidx).vf("values", values);
  return o;
}

// ------------------------------------------------------------ parameters
struct Cfg {
  ipx::Parameters p;
};

bool setParam(ipx::Parameters& p, const std::string& k, const std::string& v) {
#define IPX_PI(f) if (k == #f) { p.f = (ipxint)atol(v.c_str()); return true; }
#define IPX_PD(f) if (k == #f) { p.f = atof(v.c_str()); return true; }
  IPX_PI(dualize) IPX_PI(scale) IPX_PI(ipm_maxiter) IPX_PI(crash_basis)
  IPX_PI(rows_per_slice) IPX_PI(maxskip_updates) IPX_PI(lu_kernel)
  IPX_PI(run_crossover) IPX_PI(debug) IPX_PI(switchiter) IPX_PI(update_heuristic)
  IPX_PI(maxpasses) IPX_PI(run_centring)
  IPX_PD(ipm_feasibility_tol) IPX_PD(ipm_optimality_tol) IPX_PD(kkt_tol)
  IPX_PD(dependency_tol) IPX_PD(volume_tol) IPX_PD(lu_pivottol)
  IPX_PD(start_crossover_tol)
#undef IPX_PI
#undef IPX_PD
  return false;
}

// ipx::Parameters as IpxWrapper.cpp sets them for the HiGHS defaults
// (feasibility tol = min(primal, dual feasibility tol) = 1e-7, optimality
// 1e-8, start_crossover 1e-8, run_crossover on), silent, then the overrides.
bool makeParams(const std::string& cfg, ipx::Parameters& p) {
  p.display = 0;
  p.debug = 0;
  p.highs_logging = false;
  p.ipm_feasibility_tol = 1e-7;
  p.ipm_optimality_tol = 1e-8;
  p.start_crossover_tol = 1e-8;
  p.run_crossover = 1;
  p.time_limit = kHighsInf;
  p.ipm_maxiter = 1000;
  if (cfg == "default") return true;
  std::stringstream ss(cfg);
  std::string item;
  while (std::getline(ss, item, ',')) {
    size_t eq = item.find('=');
    if (eq == std::string::npos) return false;
    if (!setParam(p, item.substr(0, eq), item.substr(eq + 1))) {
      fprintf(stderr, "unknown ipx parameter %s\n", item.c_str());
      return false;
    }
  }
  return true;
}

Obj paramsObj(const ipx::Parameters& p, const std::string& cfg) {
  Obj o("params");
  o.s("cfg", cfg);
#define IPX_I(x) o.i(#x, (long long)p.x);
#define IPX_D(x) o.f(#x, p.x);
  IPX_I(display) IPX_D(print_interval) IPX_D(time_limit) IPX_I(analyse_basis_data)
  IPX_I(dualize) IPX_I(scale) IPX_I(ipm_maxiter) IPX_D(ipm_feasibility_tol)
  IPX_D(ipm_optimality_tol) IPX_D(ipm_drop_primal) IPX_D(ipm_drop_dual) IPX_D(kkt_tol)
  IPX_I(crash_basis) IPX_D(dependency_tol) IPX_D(volume_tol) IPX_I(rows_per_slice)
  IPX_I(maxskip_updates) IPX_I(lu_kernel) IPX_D(lu_pivottol) IPX_I(run_crossover)
  IPX_D(start_crossover_tol) IPX_D(pfeasibility_tol) IPX_D(dfeasibility_tol)
  IPX_I(debug) IPX_I(switchiter) IPX_I(stop_at_switch) IPX_I(update_heuristic)
  IPX_I(maxpasses) IPX_I(run_centring) IPX_I(max_centring_steps)
  IPX_D(centring_ratio_tolerance) IPX_D(centring_ratio_reduction)
  IPX_D(centring_alpha_scaling) IPX_I(bad_products_tolerance) IPX_I(highs_logging)
  IPX_I(timeless_log)
#undef IPX_I
#undef IPX_D
  return o;
}

// ipx::Info without the wall-clock fields (time_*).
Obj infoObj(const ipx::Info& inf) {
  Obj o;
#define IPX_I(x) o.i(#x, (long long)inf.x);
#define IPX_D(x) o.f(#x, inf.x);
  IPX_I(status) IPX_I(status_ipm) IPX_I(status_crossover) IPX_I(errflag)
  IPX_I(num_var) IPX_I(num_constr) IPX_I(num_entries) IPX_I(num_rows_solver)
  IPX_I(num_cols_solver) IPX_I(num_entries_solver) IPX_I(dualized) IPX_I(dense_cols)
  IPX_I(centring_tried) IPX_I(centring_success) IPX_I(dependent_rows)
  IPX_I(dependent_cols) IPX_I(rows_inconsistent) IPX_I(cols_inconsistent)
  IPX_I(primal_dropped) IPX_I(dual_dropped)
  IPX_D(abs_presidual) IPX_D(abs_dresidual) IPX_D(rel_presidual) IPX_D(rel_dresidual)
  IPX_D(pobjval) IPX_D(dobjval) IPX_D(rel_objgap) IPX_D(complementarity) IPX_D(normx)
  IPX_D(normy) IPX_D(normz) IPX_D(objval) IPX_D(primal_infeas) IPX_D(dual_infeas)
  IPX_I(iter) IPX_I(kktiter1) IPX_I(kktiter2) IPX_I(basis_repairs) IPX_I(updates_start)
  IPX_I(updates_ipm) IPX_I(updates_crossover)
  IPX_D(ftran_sparse) IPX_D(btran_sparse) IPX_D(mean_fill) IPX_D(max_fill)
  IPX_I(maxvol_updates) IPX_I(maxvol_skipped) IPX_I(maxvol_passes) IPX_I(tbl_nnz)
  IPX_D(tbl_max) IPX_D(frobnorm_squared) IPX_D(lambdamax) IPX_D(volume_increase)
#undef IPX_I
#undef IPX_D
  return o;
}

// ------------------------------------------------------------ model input
struct Input {
  Int num_col = 0, num_row = 0;
  double offset = 0;
  std::vector<double> obj, col_lb, col_ub, Ax, rhs;
  std::vector<Int> Ap, Ai;
  std::vector<char> constraint_type;
};

Input makeInput(const std::string& mps) {
  HighsLp lp = oracleLoadLp(mps);
  Input in;
  fillInIpxData(lp, in.num_col, in.num_row, in.offset, in.obj, in.col_lb, in.col_ub,
                in.Ap, in.Ai, in.Ax, in.rhs, in.constraint_type);
  return in;
}

Obj inputObj(const Input& in, const std::string& name) {
  Obj o("input");
  std::vector<HighsInt> ct(in.constraint_type.begin(), in.constraint_type.end());
  o.s("name", name).i("num_var", in.num_col).i("num_constr", in.num_row)
      .f("offset", in.offset).vf("obj", in.obj).vf("lb", in.col_lb).vf("ub", in.col_ub)
      .vi("Ap", in.Ap).vi("Ai", in.Ai).vf("Ax", in.Ax).vf("rhs", in.rhs)
      .vi("constr_type", ct);
  return o;
}

template <typename T>
const T* dataOrNull(const std::vector<T>& v) {
  return v.empty() ? nullptr : v.data();
}

Int loadModel(const ipx::Control& control, const Input& in, ipx::Model& model) {
  return model.Load(control, in.num_row, in.num_col, dataOrNull(in.Ap), dataOrNull(in.Ai),
                    dataOrNull(in.Ax), dataOrNull(in.rhs), dataOrNull(in.constraint_type),
                    in.offset, dataOrNull(in.obj), dataOrNull(in.col_lb),
                    dataOrNull(in.col_ub));
}

Obj modelObj(const ipx::Model& m, Int errflag) {
  Obj o("model");
  ipx::Info info;
  m.GetInfo(&info);
  std::vector<HighsInt> ct(m.constr_type_.begin(), m.constr_type_.end());
  o.i("errflag", errflag).b("dualized", m.dualized_).i("num_rows", m.num_rows_)
      .i("num_cols", m.num_cols_).i("num_dense_cols", m.num_dense_cols_)
      .i("nz_dense", m.nz_dense_).o("AI", matObj(m.AI_)).o("AIt", matObj(m.AIt_));
  vv(o, "b", m.b_);
  vv(o, "c", m.c_);
  vv(o, "lb", m.lb_);
  vv(o, "ub", m.ub_);
  o.f("norm_bounds", m.norm_bounds_).f("norm_c", m.norm_c_)
      .i("num_constr", m.num_constr_).i("num_eqconstr", m.num_eqconstr_)
      .i("num_var", m.num_var_).i("num_free_var", m.num_free_var_)
      .i("num_entries", m.num_entries_).vi("boxed_vars", m.boxed_vars_)
      .vi("constr_type", ct).f("norm_obj", m.norm_obj_).f("norm_rhs", m.norm_rhs_)
      .f("offset", m.offset_);
  vv(o, "scaled_obj", m.scaled_obj_);
  vv(o, "scaled_rhs", m.scaled_rhs_);
  vv(o, "scaled_lbuser", m.scaled_lbuser_);
  vv(o, "scaled_ubuser", m.scaled_ubuser_);
  o.o("A", matObj(m.A_)).vi("flipped_vars", m.flipped_vars_);
  vv(o, "colscale", m.colscale_);
  vv(o, "rowscale", m.rowscale_);
  o.b("filippo", m.filippoDualizationTest()).o("info", infoObj(info));
  return o;
}

// Pseudo-random vector in [lo, lo + width) from HighsRandom.
std::vector<double> rnd(HighsRandom& r, Int n, double lo, double width) {
  std::vector<double> v(n);
  for (Int i = 0; i < n; i++) v[i] = lo + width * r.fraction();
  return v;
}
Vector rndV(HighsRandom& r, Int n, double lo, double width) {
  std::vector<double> v = rnd(r, n, lo, width);
  Vector out(n);
  for (Int i = 0; i < n; i++) out[i] = v[i];
  return out;
}

// transforms: user <-> solver space of fixed vectors (model.cc 93-497).
void transforms(Out& out, const ipx::Model& model) {
  const Int n = model.num_var_, m = model.num_constr_;
  const Int ns = model.cols(), ms = model.rows();
  HighsRandom r(17);
  {
    // PresolveStartingPoint (crossover starting point).
    std::vector<double> x = rnd(r, n, -1, 2), s = rnd(r, m, -1, 2), y = rnd(r, m, -1, 2),
                        z = rnd(r, n, -1, 2);
    Vector xs(ns + ms), ys(ms), zs(ns + ms);
    model.PresolveStartingPoint(dataOrNull(x), dataOrNull(s), dataOrNull(y), dataOrNull(z),
                                xs, ys, zs);
    Obj o("presolve_starting_point");
    o.vf("x_user", x).vf("slack_user", s).vf("y_user", y).vf("z_user", z);
    vv(o, "x", xs);
    vv(o, "y", ys);
    vv(o, "z", zs);
    out.put(o);
  }
  {
    // PresolveIPMStartingPoint: nonnegative xl/xu/zl/zu, inf where the user
    // bound is infinite, as LpSolver::LoadIPMStartingPoint expects.
    std::vector<double> x = rnd(r, n, -1, 2), xl = rnd(r, n, 0, 1), xu = rnd(r, n, 0, 1),
                        s = rnd(r, m, -1, 2), y = rnd(r, m, -1, 2), zl = rnd(r, n, 0, 1),
                        zu = rnd(r, n, 0, 1);
    // User bounds: ScalePoint negates a flipped variable and swaps its
    // lower/upper parts, so its user lower bound is -scaled_ubuser_.
    for (Int j = 0; j < n; j++) {
      const bool flip = std::find(model.flipped_vars_.begin(), model.flipped_vars_.end(), j) !=
                        model.flipped_vars_.end();
      const double lbu = flip ? -model.scaled_ubuser_[j] : model.scaled_lbuser_[j];
      const double ubu = flip ? -model.scaled_lbuser_[j] : model.scaled_ubuser_[j];
      if (!std::isfinite(lbu)) { xl[j] = INFINITY; zl[j] = 0; }
      if (!std::isfinite(ubu)) { xu[j] = INFINITY; zu[j] = 0; }
    }
    // Slack and y signs by constraint type (model.cc 203-215).
    for (Int i = 0; i < m; i++) {
      const char t = model.constr_type_[i];
      if (t == '=') s[i] = 0.0;
      if (t == '<') { s[i] = std::fabs(s[i]); y[i] = -std::fabs(y[i]); }
      if (t == '>') { s[i] = -std::fabs(s[i]); y[i] = std::fabs(y[i]); }
    }
    Vector xs(ns + ms), xls(ns + ms), xus(ns + ms), ys(ms), zls(ns + ms), zus(ns + ms);
    Int err = model.PresolveIPMStartingPoint(
        dataOrNull(x), dataOrNull(xl), dataOrNull(xu), dataOrNull(s), dataOrNull(y),
        dataOrNull(zl), dataOrNull(zu), xs, xls, xus, ys, zls, zus);
    Obj o("presolve_ipm_starting_point");
    o.vf("x_user", x).vf("xl_user", xl).vf("xu_user", xu).vf("slack_user", s)
        .vf("y_user", y).vf("zl_user", zl).vf("zu_user", zu).i("errflag", err);
    vv(o, "x", xs);
    vv(o, "xl", xls);
    vv(o, "xu", xus);
    vv(o, "y", ys);
    vv(o, "zl", zls);
    vv(o, "zu", zus);
    out.put(o);
  }
  {
    // PostsolveInteriorSolution + EvaluateInteriorSolution of a solver
    // space point (positive xl/xu/zl/zu where the solver bound is finite).
    Vector x = rndV(r, ns + ms, -1, 2), xl = rndV(r, ns + ms, 0, 1),
           xu = rndV(r, ns + ms, 0, 1), y = rndV(r, ms, -1, 2),
           zl = rndV(r, ns + ms, 0, 1), zu = rndV(r, ns + ms, 0, 1);
    for (Int j = 0; j < ns + ms; j++) {
      if (!std::isfinite(model.lb(j))) { xl[j] = INFINITY; zl[j] = 0; }
      if (!std::isfinite(model.ub(j))) { xu[j] = INFINITY; zu[j] = 0; }
    }
    std::vector<double> xu_(n), xlu(n), xuu(n), su(m), yu(m), zlu(n), zuu(n);
    model.PostsolveInteriorSolution(x, xl, xu, y, zl, zu, xu_.data(), xlu.data(), xuu.data(),
                                    su.data(), yu.data(), zlu.data(), zuu.data());
    ipx::Info info;
    model.EvaluateInteriorSolution(x, xl, xu, y, zl, zu, &info);
    Obj o("postsolve_interior");
    vv(o, "x_solver", x);
    vv(o, "xl_solver", xl);
    vv(o, "xu_solver", xu);
    vv(o, "y_solver", y);
    vv(o, "zl_solver", zl);
    vv(o, "zu_solver", zu);
    o.vf("x", xu_).vf("xl", xlu).vf("xu", xuu).vf("slack", su).vf("y", yu).vf("zl", zlu)
        .vf("zu", zuu).o("info", infoObj(info));
    out.put(o);
  }
}

bool prepare(const std::string& mps, const std::string& name, const std::string& cfg,
             Out& out, ipx::Parameters& params, Input& in) {
  if (!out.ok()) return false;
  if (!makeParams(cfg, params)) return false;
  in = makeInput(mps);
  out.put(paramsObj(params, cfg));
  out.put(inputObj(in, name));
  return true;
}

// ------------------------------------------------------------ iterate
Obj iterateObj(const char* tag, const ipx::Iterate& it) {
  Obj o("iterate");
  o.s("tag", tag);
  std::vector<HighsInt> st;
  for (auto s : it.variable_state_) st.push_back((HighsInt)s);
  std::vector<HighsInt> fixed, implied, blb, bub;
  std::vector<double> scaling;
  const Int nm = (Int)it.x_.size();
  for (Int j = 0; j < nm; j++) scaling.push_back(it.ScalingFactor(j));
  vv(o, "x", it.x_);
  vv(o, "xl", it.xl_);
  vv(o, "xu", it.xu_);
  vv(o, "y", it.y_);
  vv(o, "zl", it.zl_);
  vv(o, "zu", it.zu_);
  o.vi("state", st).f("feasibility_tol", it.feasibility_tol_)
      .f("optimality_tol", it.optimality_tol_)
      .f("start_crossover_tol", it.start_crossover_tol_)
      .f("bounds_measure", it.bounds_measure_).f("costs_measure", it.costs_measure_)
      .b("postprocessed", it.postprocessed_);
  // Evaluated quantities (the getters run Evaluate()).
  o.f("pobjective", it.pobjective()).f("dobjective", it.dobjective())
      .f("pobjective_after_postproc", it.pobjective_after_postproc())
      .f("dobjective_after_postproc", it.dobjective_after_postproc())
      .f("presidual", it.presidual()).f("dresidual", it.dresidual())
      .f("complementarity", it.complementarity()).f("mu", it.mu())
      .f("mu_min", it.mu_min()).f("mu_max", it.mu_max()).f("offset", it.offset_)
      .b("feasible", it.feasible()).b("optimal", it.optimal())
      .b("term_crit_reached", it.term_crit_reached());
  vv(o, "rb", it.rb());
  vv(o, "rl", it.rl());
  vv(o, "ru", it.ru());
  vv(o, "rc", it.rc());
  o.vf("scaling", scaling);
  return o;
}

// ------------------------------------------------------------ basis
Obj factorsObj(const ipx::Basis& basis) {
  const Int m = basis.model_.rows();
  ipx::SparseMatrix L, U;
  std::vector<Int> rowperm(m), colperm(m);
  basis.GetLuFactors(&L, &U, rowperm.data(), colperm.data());
  Obj o;
  o.o("L", matObj(L)).o("U", matObj(U)).vi("rowperm", rowperm).vi("colperm", colperm);
  return o;
}

Obj basisObj(const ipx::Basis& basis) {
  Obj o;
  o.vi("basis", basis.basis_).vi("map2basis", basis.map2basis_)
      .b("fresh", basis.factorization_is_fresh_)
      .i("factorizations", basis.num_factorizations_).i("updates", basis.num_updates_)
      .i("ftran", basis.num_ftran_).i("btran", basis.num_btran_)
      .i("ftran_sparse", basis.num_ftran_sparse_).i("btran_sparse", basis.num_btran_sparse_)
      .vf("fill_factors", basis.fill_factors_)
      .f("sum_ftran_density", basis.sum_ftran_density_)
      .f("sum_btran_density", basis.sum_btran_density_)
      .f("mean_fill", basis.mean_fill()).f("max_fill", basis.max_fill())
      .f("pivottol", basis.lu_->pivottol());
  // BasicLu (lu_kernel <= 0): the BASICLU memory sizes after the Reallocate
  // loop. ForrestTomlin: the number of updates in its eta file.
  if (const ipx::BasicLu* blu = dynamic_cast<const ipx::BasicLu*>(basis.lu_.get())) {
    o.s("lu", "basiclu").i("lmem", (long long)blu->xstore_[BASICLU_MEMORYL])
        .i("umem", (long long)blu->xstore_[BASICLU_MEMORYU])
        .i("wmem", (long long)blu->xstore_[BASICLU_MEMORYW])
        .i("istore_size", (long long)blu->istore_.size())
        .i("xstore_size", (long long)blu->xstore_.size())
        .i("Li_size", (long long)blu->Li_.size()).i("Ui_size", (long long)blu->Ui_.size())
        .i("Wi_size", (long long)blu->Wi_.size())
        .f("fill_factor", blu->fill_factor_);
  } else if (const ipx::ForrestTomlin* ft =
                 dynamic_cast<const ipx::ForrestTomlin*>(basis.lu_.get())) {
    o.s("lu", "forrest_tomlin").i("ft_updates", (long long)ft->replaced_.size())
        .f("fill_factor", ft->fill_factor_);
  }
  if (basis.factorization_is_fresh_) o.o("factors", factorsObj(basis));
  return o;
}

// IPX's RandomPermute (symbolic_invert.cc:10-18, file static, so it cannot
// be called): the same code against the same libc++.
//
// The fixture pins the libc++ sequence (the generation platform, and what
// the Rust port mirrors). libstdc++ differs in both halves: its
// default_random_engine is minstd_rand0 (multiplier 16807, libc++:
// minstd_rand, 48271) and its uniform_int_distribution downscales instead
// of drawing masked independent bits. So a non-libc++ build (GCC on Linux)
// replays libc++'s algorithm here; the HiGHS library itself is untouched.
#ifndef _LIBCPP_VERSION
namespace {
// libc++ uniform_int_distribution<int>::operator()(g, {a, b}) over
// minstd_rand, through the __independent_bits_engine<minstd_rand,
// uint32_t> it builds per call. Engine range R = max - min + 1 = 2^31 - 2,
// so floor(log2 R) = 30 bits per draw; for w <= 30 bits that is one draw
// (n = 1, w0 = w), which is every IPX dimension (Int, m < 2^30).
Int libcxxUniformInt(std::minstd_rand& g, Int a, Int b) {
  const uint32_t rp = uint32_t(b) - uint32_t(a) + 1u;
  if (rp == 1) return a;
  const unsigned dt = 32;  // numeric_limits<uint32_t>::digits
  unsigned w = dt - 1;
  while (!(rp >> w)) w--;  // w = dt - countl_zero(rp) - 1
  if (rp & (0xffffffffu >> (dt - w))) ++w;
  if (w < 1 || w > 30) {
    fprintf(stderr, "highs-oracle: libcxxUniformInt range 2^%u unsupported\n", w);
    std::abort();
  }
  const uint32_t r = uint32_t(std::minstd_rand::max() - std::minstd_rand::min()) + 1u;
  const uint32_t y0 = (r >> w) << w;
  const uint32_t mask0 = 0xffffffffu >> (dt - w);
  for (;;) {
    uint32_t u;
    do {
      u = uint32_t(g() - std::minstd_rand::min());
    } while (u >= y0);
    const uint32_t s = u & mask0;
    if (s < rp) return Int(s + uint32_t(a));
  }
}
}  // namespace
#endif

std::vector<Int> randomPermute(const std::vector<Int>& basis) {
  const Int m = basis.size();
  std::vector<Int> permuted_basis(basis);
#ifdef _LIBCPP_VERSION
  std::default_random_engine re;
  std::uniform_int_distribution<Int> dist(0, m - 1);
  for (Int k = 0; k < m; k++) std::swap(permuted_basis[k], permuted_basis[dist(re)]);
#else
  std::minstd_rand re;
  for (Int k = 0; k < m; k++)
    std::swap(permuted_basis[k], permuted_basis[libcxxUniformInt(re, 0, m - 1)]);
#endif
  return permuted_basis;
}

void basisOps(Out& out, const char* tag, const ipx::Basis& basis, HighsRandom& r) {
  const Int m = basis.model_.rows();
  Obj o("basis_ops");
  o.s("tag", tag);
  Vector rhs = rndV(r, m, -1, 2);
  Vector lhs_n(m), lhs_t(m), alias(rhs);
  basis.SolveDense(rhs, lhs_n, 'N');
  basis.SolveDense(rhs, lhs_t, 'T');
  basis.SolveDense(alias, alias, 'N');  // aliased, as starting_basis.cc:76
  vv(o, "rhs", rhs);
  vv(o, "lhs_n", lhs_n);
  vv(o, "lhs_t", lhs_t);
  vv(o, "lhs_alias", alias);
  o.f("min_singular_value", basis.MinSingularValue());
  std::vector<Int> rowcounts(m), colcounts(m);
  basis.SymbolicInvert(rowcounts.data(), colcounts.data());
  o.vi("rowcounts", rowcounts).vi("colcounts", colcounts)
      .f("density_inverse", basis.DensityInverse())
      .vi("random_permute", randomPermute(basis.basis_));
  out.put(o);
}

// ------------------------------------------------------------ KKT solves
void kktSolve(Out& out, const char* type, const char* tag, ipx::KKTSolver& kkt,
              ipx::Iterate* it, ipx::Info& info, Int n, Int m, HighsRandom& r,
              double tol, const ipx::Basis* basis) {
  Obj o(type);
  o.s("tag", tag).b("with_iterate", it != nullptr).o("info_in", infoObj(info));
  kkt.Factorize(it, &info);
  o.o("info_factorize", infoObj(info)).i("basis_changes", kkt.basis_changes());
  if (it) o.o("iterate", iterateObj(tag, *it));
  if (basis) o.o("basis", basisObj(*basis));
  if (info.errflag) {
    out.put(o);
    return;
  }
  Vector a = rndV(r, n + m, -1, 2), b = rndV(r, m, -1, 2);
  Vector x(n + m), y(m);
  kkt.Solve(a, b, tol, x, y, &info);
  vv(o, "a", a);
  vv(o, "b", b);
  o.f("tol", tol);
  vv(o, "x", x);
  vv(o, "y", y);
  o.i("iter", kkt.iter()).o("info_solve", infoObj(info));
  out.put(o);
}

}  // namespace

int runIpxModel(const std::string& mps, const std::string& name, const std::string& cfg,
                const std::string& path) {
  Out out(path);
  ipx::Parameters params;
  Input in;
  if (!prepare(mps, name, cfg, out, params, in)) return 1;
  ipx::Control control;
  control.parameters(params);
  ipx::Model model;
  Int err = loadModel(control, in, model);
  out.put(modelObj(model, err));
  if (err == 0) transforms(out, model);
  return 0;
}

int runIpxLinalg(const std::string& mps, const std::string& name, const std::string& cfg,
                 const std::string& path) {
  Out out(path);
  ipx::Parameters params;
  Input in;
  if (!prepare(mps, name, cfg, out, params, in)) return 1;
  ipx::Control control;
  control.parameters(params);
  ipx::Model model;
  Int err = loadModel(control, in, model);
  out.put(modelObj(model, err));
  if (err) return 0;
  const Int m = model.rows(), n = model.cols();
  HighsRandom r(23);
  // Control::InterruptCheck polls the HiGHS task executor (Highs::run sets
  // it up; here the IPM is driven directly).
  highs::parallel::initialize_scheduler(1);

  ipx::Info info;
  model.GetInfo(&info);
  // LpSolver::InteriorPointSolve (lp_solver.cc 368-373).
  ipx::Iterate iterate(model);
  iterate.feasibility_tol(control.ipm_feasibility_tol());
  iterate.optimality_tol(control.ipm_optimality_tol());
  if (control.run_crossover()) iterate.start_crossover_tol(control.start_crossover_tol());
  out.put(iterateObj("constructed", iterate));

  ipx::IPM ipm(control);
  {
    // LpSolver::ComputeStartingPoint.
    ipx::KKTSolverDiag kkt(control, model);
    ipm.StartingPoint(&kkt, &iterate, &info);
  }
  Obj sp("ipm_starting_point");
  sp.o("info", infoObj(info));
  out.put(sp);
  if (info.status_ipm != IPX_STATUS_not_run) return 0;
  out.put(iterateObj("starting_point", iterate));
  const Int kkt_maxiter = std::min(500l, (long)(10 + m / 20));
  {
    ipx::Info scratch = info;
    ipx::KKTSolverDiag kkt(control, model);
    kkt.maxiter(kkt_maxiter);
    kktSolve(out, "kkt_diag", "starting_point_null", kkt, nullptr, scratch, n, m, r,
             0.1, nullptr);
    ipx::KKTSolverDiag kkt2(control, model);
    kkt2.maxiter(kkt_maxiter);
    kktSolve(out, "kkt_diag", "starting_point", kkt2, &iterate, scratch, n, m, r,
             control.kkt_tol() * std::sqrt(iterate.mu()), nullptr);
  }
  {
    // LpSolver::RunInitialIPM (lp_solver.cc 486-529).
    ipx::KKTSolverDiag kkt(control, model);
    if (control.switchiter() < 0) {
      kkt.maxiter(kkt_maxiter);
      ipm.maxiter(control.ipm_maxiter());
    } else {
      ipm.maxiter(std::min(control.switchiter(), control.ipm_maxiter()));
    }
    ipm.Driver(&kkt, &iterate, &info);
    switch (info.status_ipm) {
      case IPX_STATUS_optimal:
      case IPX_STATUS_no_progress:
        info.status_ipm = IPX_STATUS_not_run;
        break;
      case IPX_STATUS_failed:
        info.status_ipm = IPX_STATUS_not_run;
        info.errflag = 0;
        break;
      case IPX_STATUS_iter_limit:
        if (info.iter < control.ipm_maxiter()) info.status_ipm = IPX_STATUS_not_run;
    }
  }
  Obj ii("ipm_initial");
  ii.o("info", infoObj(info));
  out.put(ii);
  out.put(iterateObj("initial_ipm", iterate));
  {
    ipx::Info scratch = info;
    ipx::KKTSolverDiag kkt(control, model);
    kkt.maxiter(kkt_maxiter);
    kktSolve(out, "kkt_diag", "initial_ipm", kkt, &iterate, scratch, n, m, r,
             control.kkt_tol() * std::sqrt(iterate.mu()), nullptr);
  }

  // LpSolver::BuildStartingBasis.
  ipx::Basis basis(control, model);
  {
    Obj o("starting_basis");
    o.o("info_in", infoObj(info));
    ipx::StartingBasis(&iterate, &basis, &info);
    o.o("info", infoObj(info)).o("basis", basisObj(basis))
        .o("iterate", iterateObj("starting_basis", iterate));
    out.put(o);
  }
  if (info.errflag) return 0;
  if (model.dualized()) {
    std::swap(info.dependent_rows, info.dependent_cols);
    std::swap(info.rows_inconsistent, info.cols_inconsistent);
  }
  basisOps(out, "starting_basis", basis, r);

  // KKTSolverBasis on the starting basis, then a few main IPM iterations.
  {
    ipx::KKTSolverBasis kkt(control, basis);
    kktSolve(out, "kkt_basis", "starting_basis", kkt, &iterate, info, n, m, r,
             control.kkt_tol() * std::sqrt(iterate.mu()), &basis);
  }
  if (info.errflag) return 0;
  {
    ipx::KKTSolverBasis kkt(control, basis);
    ipm.maxiter(std::min(info.iter + 4, control.ipm_maxiter()));
    ipm.Driver(&kkt, &iterate, &info);
  }
  Obj mi("ipm_main");
  mi.o("info", infoObj(info));
  out.put(mi);
  out.put(iterateObj("main_ipm", iterate));
  std::vector<int> status(n + m);
  for (Int j = 0; j < n + m; j++) status[j] = basis.StatusOf(j);

  // A fresh Basis loaded with the statuses (both sides start from the same
  // factorization), then KKTSolverBasis on the late iterate.
  ipx::Basis basis2(control, model);
  {
    Obj o("basis_load");
    std::vector<HighsInt> st(status.begin(), status.end());
    o.vi("status", st);
    Int lerr = loadBasis(basis2, status.data());
    o.i("errflag", lerr).o("basis", basisObj(basis2));
    out.put(o);
    if (lerr) return 0;
  }
  basisOps(out, "basis_load", basis2, r);
  {
    ipx::KKTSolverBasis kkt(control, basis2);
    kktSolve(out, "kkt_basis", "basis_load", kkt, &iterate, info, n, m, r,
             control.kkt_tol() * std::sqrt(iterate.mu()), &basis2);
  }

  // Postprocessing and postsolve (lp_solver.cc 377-378, 568-584).
  {
    Obj o("postsolve");
    // LpSolver::InteriorPointSolve: Postprocess + EvaluatePostsolved right
    // after the IPM, then GetInteriorSolution postsolves the postprocessed
    // iterate (implied variables have finite xl/xu again).
    iterate.Postprocess();
    ipx::Info pinfo = info;
    iterate.EvaluatePostsolved(&pinfo);
    o.o("info", infoObj(pinfo)).o("iterate", iterateObj("postprocessed", iterate));
    const Int nu = model.num_var_, mu = model.num_constr_;
    std::vector<double> xu(nu), xlu(nu), xuu(nu), su(mu), yu(mu), zlu(nu), zuu(nu);
    model.PostsolveInteriorSolution(iterate.x(), iterate.xl(), iterate.xu(), iterate.y(),
                                    iterate.zl(), iterate.zu(), xu.data(), xlu.data(),
                                    xuu.data(), su.data(), yu.data(), zlu.data(), zuu.data());
    o.vf("x", xu).vf("xl", xlu).vf("xu", xuu).vf("slack", su).vf("y", yu).vf("zl", zlu)
        .vf("zu", zuu);
    Vector xc(n + m), yc(m), zc(n + m);
    iterate.DropToComplementarity(xc, yc, zc);
    vv(o, "x_crossover", xc);
    vv(o, "y_crossover", yc);
    vv(o, "z_crossover", zc);
    // Basic statuses for this basis: nonbasic at the bound x sits on,
    // superbasic otherwise; a slack column (j >= n) is never superbasic
    // (asserted in Model::DualizeBackBasis): nonbasic at its finite bound.
    std::vector<Int> basic_status(n + m);
    for (Int j = 0; j < n + m; j++)
      basic_status[j] = basis2.IsBasic(j)                     ? IPX_basic
                        : xc[j] == model.lb(j)                ? IPX_nonbasic_lb
                        : xc[j] == model.ub(j)                ? IPX_nonbasic_ub
                        : j < n                               ? IPX_superbasic
                        : std::isfinite(model.lb(j))          ? IPX_nonbasic_lb
                        : std::isfinite(model.ub(j))          ? IPX_nonbasic_ub
                                                              : IPX_nonbasic;
    Vector xb(xc), yb(yc), zb(zc);
    basis2.ComputeBasicSolution(xb, yb, zb);
    vv(o, "x_basic", xb);
    vv(o, "y_basic", yb);
    vv(o, "z_basic", zb);
    o.vi("basic_status", basic_status);
    std::vector<double> xbu(nu), sbu(mu), ybu(mu), zbu(nu);
    model.PostsolveBasicSolution(xb, yb, zb, basic_status, xbu.data(), sbu.data(),
                                 ybu.data(), zbu.data());
    o.vf("xb_user", xbu).vf("slackb_user", sbu).vf("yb_user", ybu).vf("zb_user", zbu);
    ipx::Info binfo;
    model.EvaluateBasicSolution(xb, yb, zb, basic_status, &binfo);
    o.o("info_basic", infoObj(binfo));
    std::vector<Int> cbasis(mu), vbasis(nu);
    model.PostsolveBasis(basic_status, cbasis.data(), vbasis.data());
    o.vi("cbasis", cbasis).vi("vbasis", vbasis);
    out.put(o);
  }
  return 0;
}

// ipx-trace <mps> <name> <cfg> <out>: oracle (b), the per-iteration IPX trace
// of a full Highs::run with solver=ipx (presolve off, serial). The records
// are written by the second patch, highs-oracle-ipx-hooks.patch
// (HighsOracleIpxDump.h, FORMAT.md "ipx/trace"), into a temporary file
// that is then framed by a trace_header and a result record. <cfg> is a
// comma separated list of HiGHS option overrides (key=value, or
// "default").
int runIpxTrace(const std::string& mps, const std::string& name, const std::string& cfg,
                const std::string& path) {
  Highs highs;
  oracleBaseOptions(highs, 0);
  highs.setOptionValue("solver", "ipx");
  if (cfg != "default") {
    std::stringstream ss(cfg);
    std::string kv;
    while (std::getline(ss, kv, ',')) {
      const size_t eq = kv.find('=');
      if (eq == std::string::npos ||
          highs.setOptionValue(kv.substr(0, eq), kv.substr(eq + 1)) != HighsStatus::kOk) {
        fprintf(stderr, "highs-oracle: bad option %s\n", kv.c_str());
        return 64;
      }
    }
  }
  if (highs.readModel(mps) == HighsStatus::kError) {
    fprintf(stderr, "highs-oracle: cannot read %s\n", mps.c_str());
    return 2;
  }
  highs.clearIntegrality();
  const std::string tmp = path + ".hook.tmp";
  std::remove(tmp.c_str());
  setenv("HIGHS_ORACLE_IPX_DUMP", tmp.c_str(), 1);
  HighsStatus run_status = highs.run();
  unsetenv("HIGHS_ORACLE_IPX_DUMP");
  std::FILE* hook = ipx::highsOracleIpxFile();
  if (hook == nullptr) {
    fprintf(stderr, "highs-oracle: IPX hook did not run (is highs-oracle-ipx-hooks.patch applied?)\n");
    return 3;
  }
  std::fflush(hook);
  std::string hook_lines;
  {
    std::ifstream in(tmp, std::ios::binary);
    std::stringstream buf;
    buf << in.rdbuf();
    hook_lines = buf.str();
  }
  std::remove(tmp.c_str());

  const HighsInfo& info = highs.getInfo();
  const HighsSolution& sol = highs.getSolution();
  const HighsBasis& basis = highs.getBasis();
  std::vector<HighsInt> col_status, row_status;
  for (auto s : basis.col_status) col_status.push_back((HighsInt)s);
  for (auto s : basis.row_status) row_status.push_back((HighsInt)s);
  Obj header("trace_header");
  header.s("name", name).s("cfg", cfg).s("solver", "ipx")
      .s("highs_version", highsVersion()).s("highs_githash", highsGithash());
  Obj result("result");
  result.i("run_status", (int)run_status)
      .i("model_status", (int)highs.getModelStatus())
      .s("model_status_string", highs.modelStatusToString(highs.getModelStatus()))
      .f("objective_function_value", info.objective_function_value)
      .i("ipm_iteration_count", info.ipm_iteration_count)
      .i("crossover_iteration_count", info.crossover_iteration_count)
      .i("simplex_iteration_count", info.simplex_iteration_count)
      .i("primal_solution_status", info.primal_solution_status)
      .i("dual_solution_status", info.dual_solution_status)
      .b("basis_valid", basis.valid).vi("col_status", col_status)
      .vi("row_status", row_status).vf("col_value", sol.col_value)
      .vf("col_dual", sol.col_dual).vf("row_value", sol.row_value)
      .vf("row_dual", sol.row_dual);
  std::FILE* f = std::fopen(path.c_str(), "wb");
  if (f == nullptr) {
    fprintf(stderr, "highs-oracle: cannot open %s\n", path.c_str());
    return 2;
  }
  const std::string h = header.str() + "\n", r = result.str() + "\n";
  std::fwrite(h.data(), 1, h.size(), f);
  std::fwrite(hook_lines.data(), 1, hook_lines.size(), f);
  std::fwrite(r.data(), 1, r.size(), f);
  std::fclose(f);
  return 0;
}

// ipx <mps> <name> <out>: oracle (c), end-to-end IPX through Highs::run
// (solver=ipx, presolve off, serial) in several configurations on one LP.
// Records (FORMAT.md "ipx/run"):
//   ipx_header  name, versions, hipo_available (the useHipo probe) and the LP HiGHS read (lp, unscaled).
//   per cfg:    run_header {cfg, options overridden}, the oracle (b) hook
//               records of that run (call counter restarted at 0), and
//               run_result: run/model status, the whole HighsInfo (LP part),
//               HighsSolution and HighsBasis.
// Configurations (comma separated HiGHS option overrides on top of
// solver=ipx): xon (run_crossover=on), xoff (run_crossover=off), ipm
// (solver=ipm: the routing of "ipm" when HiPO is unavailable), dual
// (ipx_dualize_strategy=1, crossover on).
namespace {

Obj ipxLpObj(const HighsLp& lp) {
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

Obj highsInfoObj(const HighsInfo& info) {
  Obj o;
#define HI_I(x) o.i(#x, (long long)info.x);
#define HI_D(x) o.f(#x, info.x);
  HI_I(valid) HI_I(simplex_iteration_count) HI_I(ipm_iteration_count)
  HI_I(crossover_iteration_count) HI_I(pdlp_iteration_count) HI_I(qp_iteration_count)
  HI_I(primal_solution_status) HI_I(dual_solution_status) HI_I(basis_validity)
  HI_D(objective_function_value)
  HI_I(num_primal_infeasibilities) HI_D(max_primal_infeasibility)
  HI_D(sum_primal_infeasibilities) HI_I(num_dual_infeasibilities)
  HI_D(max_dual_infeasibility) HI_D(sum_dual_infeasibilities)
  HI_I(num_relative_primal_infeasibilities) HI_D(max_relative_primal_infeasibility)
  HI_I(num_relative_dual_infeasibilities) HI_D(max_relative_dual_infeasibility)
  HI_I(num_primal_residual_errors) HI_D(max_primal_residual_error)
  HI_I(num_dual_residual_errors) HI_D(max_dual_residual_error)
  HI_I(num_relative_primal_residual_errors) HI_D(max_relative_primal_residual_error)
  HI_I(num_relative_dual_residual_errors) HI_D(max_relative_dual_residual_error)
  HI_I(num_complementarity_violations) HI_D(max_complementarity_violation)
  HI_D(primal_dual_objective_error)
#undef HI_I
#undef HI_D
  return o;
}

bool applyOverrides(Highs& highs, const std::string& cfg) {
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

}  // namespace

int runIpx(const std::string& mps, const std::string& name, const std::string& path) {
  static const char* const kCfgs[][2] = {
      {"xon", "run_crossover=on"},
      {"xoff", "run_crossover=off"},
      {"ipm", "solver=ipm,run_crossover=on"},
      {"dual", "ipx_dualize_strategy=1,run_crossover=on"}};
  HighsLp lp;
  {
    Highs reader;
    oracleBaseOptions(reader, 0);
    if (reader.readModel(mps) == HighsStatus::kError) {
      fprintf(stderr, "highs-oracle: cannot read %s\n", mps.c_str());
      return 2;
    }
    reader.clearIntegrality();
    lp = reader.getLp();
  }
  Out out(path);
  if (!out.ok()) return 2;
  Obj header("ipx_header");
  header.s("name", name).s("highs_version", highsVersion())
      .s("highs_githash", highsGithash())
      .b("hipo_available", HighsExternalApi::isAvailable<HighsExtras::hipo>())
      .o("lp", ipxLpObj(lp));
  out.put(header);
  const std::string tmp = path + ".hook.tmp";
  std::remove(tmp.c_str());
  setenv("HIGHS_ORACLE_IPX_DUMP", tmp.c_str(), 1);
  for (const auto& c : kCfgs) {
    Highs highs;
    oracleBaseOptions(highs, 0);
    highs.setOptionValue("solver", "ipx");
    if (!applyOverrides(highs, c[1])) return 64;
    Obj rh("run_header");
    rh.s("cfg", c[0]).s("overrides", c[1]);
    out.put(rh);
    out.flush();
    if (highs.passModel(lp) != HighsStatus::kOk) {
      fprintf(stderr, "highs-oracle: passModel failed for %s\n", name.c_str());
      return 2;
    }
    // Restart the hook's call counter and empty its file for this run.
    ipx::highsOracleIpxCall() = 0;
    std::FILE* hook = ipx::highsOracleIpxFile();
    if (hook == nullptr) {
      fprintf(stderr, "highs-oracle: IPX hook file not open (is highs-oracle-ipx-hooks.patch applied?)\n");
      return 3;
    }
    std::fflush(hook);
    if (ftruncate(fileno(hook), 0) != 0) return 3;
    std::rewind(hook);
    const HighsStatus run_status = highs.run();
    std::fflush(hook);
    {
      std::ifstream in(tmp, std::ios::binary);
      std::stringstream buf;
      buf << in.rdbuf();
      out.raw_lines(buf.str());
    }
    const HighsSolution& sol = highs.getSolution();
    const HighsBasis& basis = highs.getBasis();
    std::vector<HighsInt> col_status, row_status;
    for (auto s : basis.col_status) col_status.push_back((HighsInt)s);
    for (auto s : basis.row_status) row_status.push_back((HighsInt)s);
    Obj r("run_result");
    r.s("cfg", c[0]).i("run_status", (int)run_status)
        .i("model_status", (int)highs.getModelStatus())
        .o("info", highsInfoObj(highs.getInfo()))
        .b("value_valid", sol.value_valid).b("dual_valid", sol.dual_valid)
        .vf("col_value", sol.col_value).vf("col_dual", sol.col_dual)
        .vf("row_value", sol.row_value).vf("row_dual", sol.row_dual)
        .b("basis_valid", basis.valid).b("basis_alien", basis.alien)
        .vi("col_status", col_status).vi("row_status", row_status);
    out.put(r);
  }
  unsetenv("HIGHS_ORACLE_IPX_DUMP");
  std::remove(tmp.c_str());
  return 0;
}
