// HFactor fixtures: build (INVERT) on several bases of the *scaled* LP that
// HEkk factorises, FTRAN/BTRAN on dense / sparse / hyper-sparse RHS, update
// sequences under FT / PF / MPF / APF with post-update solves, and the
// HFactorRefactor path (refactor_info_.use = true).
#include <cmath>

#include "lp_data/HConst.h"
#include "oracle.h"
#include "util/HFactor.h"
#include "util/HighsRandom.h"

namespace {

const double kSolveDensities[3] = {1.0, 0.1, 0.01};

HVector makeVec(HighsInt size) {
  HVector v;
  v.setup(size);
  v.clear();
  return v;
}

void setupFactor(HFactor& f, const HighsSparseMatrix& a,
                 std::vector<HighsInt>& basic_index, HighsInt update_method) {
  // Mirrors HSimplexNla::setup (num_basic = num_row, default thresholds,
  // use_original_HFactor_logic = true), plus an explicit update method.
  f.setupGeneral(a.num_col_, a.num_row_, a.num_row_, a.start_.data(),
                 a.index_.data(), a.value_.data(), basic_index.data(),
                 kDefaultPivotThreshold, kDefaultPivotTolerance,
                 kHighsDebugLevelMin, NULL, true, update_method);
}

Obj invertObj(const HFactor& f) {
  InvertibleRepresentation inv = f.getInvert();
  Obj o;
  o.vi("l_pivot_index", inv.l_pivot_index).vi("l_pivot_lookup", inv.l_pivot_lookup)
      .vi("l_start", inv.l_start).vi("l_index", inv.l_index).vf("l_value", inv.l_value)
      .vi("lr_start", inv.lr_start).vi("lr_index", inv.lr_index)
      .vf("lr_value", inv.lr_value).vi("u_pivot_lookup", inv.u_pivot_lookup)
      .vi("u_pivot_index", inv.u_pivot_index).vf("u_pivot_value", inv.u_pivot_value)
      .vi("u_start", inv.u_start).vi("u_last_p", inv.u_last_p)
      .vi("u_index", inv.u_index).vf("u_value", inv.u_value)
      .vi("ur_start", inv.ur_start).vi("ur_lastp", inv.ur_lastp)
      .vi("ur_space", inv.ur_space).vi("ur_index", inv.ur_index)
      .vf("ur_value", inv.ur_value).vi("pf_start", inv.pf_start)
      .vi("pf_index", inv.pf_index).vf("pf_value", inv.pf_value)
      .vi("pf_pivot_index", inv.pf_pivot_index)
      .vf("pf_pivot_value", inv.pf_pivot_value);
  return o;
}

Obj buildObj(const HFactor& f, HighsInt rank_deficiency,
             const std::vector<HighsInt>& basic_after) {
  Obj o;
  o.i("rank_deficiency", rank_deficiency).vi("basic_index_after", basic_after)
      .vi("row_with_no_pivot", f.row_with_no_pivot)
      .vi("col_with_no_pivot", f.col_with_no_pivot)
      .vi("var_with_no_pivot", f.var_with_no_pivot)
      .f("build_synthetic_tick", f.build_synthetic_tick)
      .i("basis_matrix_num_el", f.basis_matrix_num_el)
      .i("invert_num_el", f.invert_num_el).i("kernel_dim", f.kernel_dim)
      .i("kernel_num_el", f.kernel_num_el)
      .vi("refactor_pivot_row", f.refactor_info_.pivot_row)
      .vi("refactor_pivot_var", f.refactor_info_.pivot_var)
      .vi8("refactor_pivot_type", f.refactor_info_.pivot_type)
      .f("refactor_build_synthetic_tick", f.refactor_info_.build_synthetic_tick)
      .o("invert", invertObj(f));
  return o;
}

// The RHS set used for every solve batch: three unit vectors, a structural
// column (as col_aq is formed), a sparse random and a dense random vector.
std::vector<std::pair<std::string, HVector>> solveRhs(const HighsSparseMatrix& a,
                                                      HighsUInt seed) {
  const HighsInt m = a.num_row_;
  HighsRandom rng(seed);
  std::vector<std::pair<std::string, HVector>> out;
  const HighsInt units[3] = {0, m / 2, m - 1};
  for (HighsInt r : units) {
    HVector v = makeVec(m);
    v.index[v.count++] = r;
    v.array[r] = 1.0;
    out.emplace_back("unit", v);
  }
  {
    HVector v = makeVec(m);
    a.collectAj(v, a.num_col_ / 2, 1.0);
    out.emplace_back("column", v);
  }
  {
    HVector v = makeVec(m);
    for (HighsInt r = 1; r < m; r += 7) {
      v.index[v.count++] = r;
      v.array[r] = 2 * rng.fraction() - 1;
    }
    out.emplace_back("sparse", v);
  }
  {
    HVector v = makeVec(m);
    for (HighsInt r = 0; r < m; r++) {
      v.index[v.count++] = r;
      v.array[r] = 2 * rng.fraction() - 1;
    }
    out.emplace_back("dense", v);
  }
  return out;
}

// full: every RHS x every expected density, FTRAN and BTRAN, plus the
// std::vector<double> overloads (count = -1 path) for the dense RHS.
// !full (after each update): the middle unit RHS at density 0.01 and the
// dense RHS at density 1.0, FTRAN and BTRAN.
std::string solveBatch(HFactor& f, const HighsSparseMatrix& a, HighsUInt seed,
                       bool full) {
  std::string arr = "[";
  bool first = true;
  HighsInt unit_k = 0;
  for (auto& named : solveRhs(a, seed)) {
    const bool is_unit = named.first == "unit";
    if (is_unit) unit_k++;
    if (!full && !(named.first == "dense" || (is_unit && unit_k == 2))) continue;
    for (double ed : kSolveDensities) {
      if (!full && ed != (is_unit ? 0.01 : 1.0)) continue;
      for (int dir = 0; dir < 2; dir++) {
        HVector v = named.second;
        if (dir == 0)
          f.ftranCall(v, ed);
        else
          f.btranCall(v, ed);
        Obj o;
        o.s("op", dir == 0 ? "ftran" : "btran").s("rhs_kind", named.first)
            .f("expected_density", ed).hv("in", named.second).hv("out", v);
        if (!first) arr += ',';
        first = false;
        arr += o.str();
      }
    }
    if (full && named.first == "dense") {
      for (int dir = 0; dir < 2; dir++) {
        std::vector<double> x = named.second.array;
        if (dir == 0)
          f.ftranCall(x);
        else
          f.btranCall(x);
        Obj o;
        o.s("op", dir == 0 ? "ftran_vec" : "btran_vec").s("rhs_kind", named.first)
            .vf("in", named.second.array).vf("out", x);
        arr += ',';
        arr += o.str();
      }
    }
  }
  arr += ']';
  return arr;
}

const char* methodName(HighsInt m) {
  switch (m) {
    case kUpdateMethodFt: return "FT";
    case kUpdateMethodPf: return "PF";
    case kUpdateMethodMpf: return "MPF";
    case kUpdateMethodApf: return "APF";
  }
  return "?";
}

// A deterministic simplex-like basis change sequence. Entering q: random
// nonbasic variable; aq = B^{-1} a_q formed exactly as HEkkDual::updateFtran
// does (clear, packFlag, collectAj, ftran); leaving row p = first argmax
// |aq|, rejected below 1e-7; ep = B^{-T} e_p as in HEkkDual::chooseRow.
void updateSequence(Out& out, const std::string& name, const std::string& basis_kind,
                    const HighsSparseMatrix& a, std::vector<HighsInt> basic_index,
                    HighsInt method, HighsInt num_updates) {
  const HighsInt m = a.num_row_, n_tot = a.num_col_ + a.num_row_;
  HFactor f;
  setupFactor(f, a, basic_index, method);
  std::vector<HighsInt> basic_in = basic_index;
  HighsInt rank_deficiency = f.build();
  out.put(Obj("update_sequence_start").s("lp", name).s("basis", basis_kind)
              .s("method", methodName(method)).i("update_method", method)
              .vi("basic_index", basic_in).i("rank_deficiency", rank_deficiency)
              .vi("basic_index_after", basic_index));
  if (rank_deficiency) return;
  std::vector<int8_t> is_basic(n_tot, 0);
  for (HighsInt v : basic_index) is_basic[v] = 1;
  HighsRandom rng(1000 + method);
  const double aq_density = 0.1, ep_density = 0.1;
  for (HighsInt k = 0; k < num_updates; k++) {
    HighsInt q = -1, p = -1;
    HVector aq = makeVec(m);
    for (HighsInt tries = 0; tries < 200; tries++) {
      HighsInt cand = rng.integer(n_tot);
      if (is_basic[cand]) continue;
      aq.clear();
      aq.packFlag = true;
      a.collectAj(aq, cand, 1.0);
      f.ftranCall(aq, aq_density);
      double best = 0;
      HighsInt best_i = -1;
      for (HighsInt i = 0; i < aq.count; i++) {
        HighsInt r = aq.index[i];
        double v = std::fabs(aq.array[r]);
        if (v > best || (v == best && best_i >= 0 && r < best_i)) {
          best = v;
          best_i = r;
        }
      }
      if (best_i >= 0 && best >= 1e-7) {
        q = cand;
        p = best_i;
        break;
      }
    }
    if (q < 0) break;
    HVector ep = makeVec(m);
    ep.clear();
    ep.count = 1;
    ep.index[0] = p;
    ep.array[p] = 1;
    ep.packFlag = true;
    f.btranCall(ep, ep_density);
    HighsInt hint = 0;
    HighsInt row = p;
    Obj o("update_step");
    o.s("lp", name).s("basis", basis_kind).s("method", methodName(method))
        .i("k", k).i("variable_in", q).i("row_out", p)
        .i("variable_out", basic_index[p]).f("aq_density", aq_density)
        .f("ep_density", ep_density).hv("aq", aq, true).hv("ep", ep, true);
    f.update(&aq, &ep, &row, &hint);
    is_basic[basic_index[p]] = 0;
    is_basic[q] = 1;
    basic_index[p] = q;
    o.i("hint", hint).i("row_after", row).hv("aq_after", aq, true)
        .hv("ep_after", ep, true)
        .raw("solves", solveBatch(f, a, 77 + (HighsUInt)k, false));
    out.put(o);
    if (hint) break;
  }
  out.put(Obj("update_sequence_end").s("lp", name).s("basis", basis_kind)
              .s("method", methodName(method)).vi("basic_index", basic_index)
              .o("invert", invertObj(f)));
}

void buildCase(Out& out, const std::string& name, const std::string& basis_kind,
               const HighsSparseMatrix& a, const std::vector<HighsInt>& basis) {
  std::vector<HighsInt> basic_index = basis;
  HFactor f;
  setupFactor(f, a, basic_index, kUpdateMethodFt);
  HighsInt rank_deficiency = f.build();
  out.put(Obj("build").s("lp", name).s("basis", basis_kind)
              .vi("basic_index", basis)
              .o("result", buildObj(f, rank_deficiency, basic_index))
              .raw("solves", solveBatch(f, a, 5, true)));
  if (rank_deficiency) return;
  // HFactorRefactor: a fresh factor given the recorded pivot sequence with
  // refactor_info_.use = true takes HFactor::rebuild instead of INVERT.
  std::vector<HighsInt> basic_index2 = basis;
  HFactor g;
  setupFactor(g, a, basic_index2, kUpdateMethodFt);
  g.refactor_info_ = f.refactor_info_;
  g.refactor_info_.use = true;
  HighsInt rd2 = g.build();
  out.put(Obj("refactor").s("lp", name).s("basis", basis_kind)
              .vi("basic_index", basis).vi("refactor_pivot_row", f.refactor_info_.pivot_row)
              .vi("refactor_pivot_var", f.refactor_info_.pivot_var)
              .vi8("refactor_pivot_type", f.refactor_info_.pivot_type)
              .f("refactor_build_synthetic_tick", f.refactor_info_.build_synthetic_tick)
              .o("result", buildObj(g, rd2, basic_index2))
              .raw("solves", solveBatch(g, a, 5, true)));
}

}  // namespace

int runFactor(const std::string& mps, const std::string& name,
              const std::string& path) {
  EkkSnapshot opt = oracleSnapshot(mps, -1);
  if (!opt.valid) {
    fprintf(stderr, "highs-oracle factor: no HEkk snapshot for %s\n", mps.c_str());
    return 1;
  }
  const HighsSparseMatrix& a = opt.scaled_lp.a_matrix_;
  const HighsInt num_col = a.num_col_, num_row = a.num_row_;
  EkkSnapshot mid = oracleSnapshot(mps, opt.iteration_count / 2);

  Out out(path);
  out.put(Obj("factor_lp").s("name", name).i("num_col", num_col)
              .i("num_row", num_row).i("optimal_iteration_count", opt.iteration_count)
              .i("mid_iteration_count", mid.iteration_count)
              .f("pivot_threshold", kDefaultPivotThreshold)
              .f("pivot_tolerance", kDefaultPivotTolerance).mat("a", a));

  std::vector<HighsInt> slack(num_row);
  for (HighsInt i = 0; i < num_row; i++) slack[i] = num_col + i;

  // Rank-deficient #1: the optimal basis with every third basic variable
  // replaced by a copy of the first structural basic variable.
  std::vector<HighsInt> dup = opt.basic_index;
  HighsInt first_struct = -1;
  for (HighsInt v : dup)
    if (v < num_col) {
      first_struct = v;
      break;
    }
  if (first_struct < 0) first_struct = 0;
  for (HighsInt i = 1; i < num_row; i += 3) dup[i] = first_struct;
  // Rank-deficient #2: logicals except the first rows taken by the first
  // structural columns, with two rows sharing the same structural.
  std::vector<HighsInt> mixed = slack;
  for (HighsInt i = 0; i < num_row && i < num_col; i += 2) mixed[i] = i / 2;
  if (num_row > 3) mixed[3] = mixed[0];

  buildCase(out, name, "slack", a, slack);
  buildCase(out, name, "mid", a, mid.basic_index);
  buildCase(out, name, "optimal", a, opt.basic_index);
  buildCase(out, name, "rank_deficient_dup", a, dup);
  buildCase(out, name, "rank_deficient_mixed", a, mixed);

  const HighsInt methods[4] = {kUpdateMethodFt, kUpdateMethodPf, kUpdateMethodMpf,
                               kUpdateMethodApf};
  const HighsInt num_updates = 24;
  for (HighsInt method : methods) {
    updateSequence(out, name, "slack", a, slack, method, num_updates);
    updateSequence(out, name, "mid", a, mid.basic_index, method, num_updates);
  }
  return 0;
}
