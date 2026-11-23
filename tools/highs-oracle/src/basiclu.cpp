// BASICLU fixtures: the LU kernel of
// IPX (highs/ipm/basiclu, (C) 2016-2018 ERGO-Code, distributed with HiGHS
// under the MIT license) driven through its object API on bases of the
// *scaled* LP HEkk factorises (the same bases as factor.cpp):
//
//   factorize       basiclu_obj_factorize; xstore scalar block, factors from
//                   basiclu_obj_get_factors, dense solves 'N'/'T' and sparse
//                   solves 'N'/'T'
//   update_step     one Forrest-Tomlin step: solve_for_update 'N' of the
//                   entering column and 'T' of the leaving position, then
//                   basiclu_obj_update(xtbl); xstore and solves after it
//   refactorize     basiclu_obj_factorize of the final updated basis on the
//                   same (used) store
//
// Every record carries its inputs, so the Rust replay needs no MPS reader.
#include <cmath>
#include <cstring>

#include "ipm/basiclu/basiclu.h"
#include "oracle.h"
#include "util/HighsRandom.h"

namespace {

// Column-wise basis matrix B = [A I](:, basic_index): structural j is column
// j of A, logical n+i is the unit column e_i (HiGHS' logical convention).
struct BasisMatrix {
  std::vector<lu_int> begin, end, index;
  std::vector<double> value;
};

BasisMatrix basisMatrix(const HighsSparseMatrix& a,
                        const std::vector<HighsInt>& basic_index) {
  BasisMatrix b;
  const HighsInt n = a.num_col_;
  for (HighsInt var : basic_index) {
    b.begin.push_back((lu_int)b.index.size());
    if (var < n) {
      for (HighsInt p = a.start_[var]; p < a.start_[var + 1]; p++) {
        b.index.push_back(a.index_[p]);
        b.value.push_back(a.value_[p]);
      }
    } else {
      b.index.push_back(var - n);
      b.value.push_back(1.0);
    }
    b.end.push_back((lu_int)b.index.size());
  }
  return b;
}

// basiclu_obj_initialize mallocs istore/xstore, so the slots that
// lu_initialize does not write (unused scalar slots, the half of work0 that
// its `memset(work0, 0, m * sizeof(lu_int))` misses, the update/solve
// workspaces) are indeterminate. IPX's BasicLu (basiclu_wrapper.cc) keeps
// them in zero-initialised std::vectors, and the Rust port in zeroed Vecs.
// So the oracle zeroes the fresh object and runs basiclu_initialize again:
// this is basiclu_obj_initialize with calloc in place of malloc.
void objInitZeroed(struct basiclu_object* obj, lu_int m) {
  basiclu_obj_initialize(obj, m);
  memset(obj->istore, 0, (BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * m) * sizeof(lu_int));
  memset(obj->xstore, 0, (BASICLU_SIZE_XSTORE_1 + BASICLU_SIZE_XSTORE_M * m) * sizeof(double));
  memset(obj->Li, 0, m * sizeof(lu_int));
  memset(obj->Lx, 0, m * sizeof(double));
  memset(obj->Ui, 0, m * sizeof(lu_int));
  memset(obj->Ux, 0, m * sizeof(double));
  memset(obj->Wi, 0, m * sizeof(lu_int));
  memset(obj->Wx, 0, m * sizeof(double));
  memset(obj->ilhs, 0, m * sizeof(lu_int));
  basiclu_initialize(m, obj->istore, obj->xstore);
  obj->xstore[BASICLU_MEMORYL] = m;
  obj->xstore[BASICLU_MEMORYU] = m;
  obj->xstore[BASICLU_MEMORYW] = m;
}

// xstore[0..512) (the scalar block) as the positions with a non-zero bit
// pattern and their values.
Obj xheadObj(const double* xstore) {
  std::vector<HighsInt> idx;
  std::vector<double> val;
  for (HighsInt k = 0; k < 512; k++) {
    uint64_t bits;
    memcpy(&bits, &xstore[k], 8);
    if (bits) {
      idx.push_back(k);
      val.push_back(xstore[k]);
    }
  }
  Obj o;
  o.vi("i", idx).vf("v", val);
  return o;
}

Obj factorsObj(struct basiclu_object* obj, lu_int m) {
  const lu_int lnz = (lu_int)obj->xstore[BASICLU_LNZ];
  const lu_int unz = (lu_int)obj->xstore[BASICLU_UNZ];
  std::vector<lu_int> rowperm(m), colperm(m), lcolptr(m + 1), lrowidx(m + lnz),
      ucolptr(m + 1), urowidx(m + unz);
  std::vector<double> lvalue(m + lnz), uvalue(m + unz);
  lu_int status = basiclu_obj_get_factors(obj, rowperm.data(), colperm.data(),
                                          lcolptr.data(), lrowidx.data(), lvalue.data(),
                                          ucolptr.data(), urowidx.data(), uvalue.data());
  Obj o;
  o.i("status", status).vi("rowperm", rowperm).vi("colperm", colperm)
      .vi("lcolptr", lcolptr).vi("lrowidx", lrowidx).vf("lvalue", lvalue)
      .vi("ucolptr", ucolptr).vi("urowidx", urowidx).vf("uvalue", uvalue)
      .o("xhead", xheadObj(obj->xstore));
  return o;
}

// A length-m vector as the positions with a non-zero bit pattern and their
// values (so -0.0 survives; every other entry is +0.0).
Obj nzObj(const double* x, lu_int m) {
  std::vector<HighsInt> idx;
  std::vector<double> val;
  for (lu_int k = 0; k < m; k++) {
    uint64_t bits;
    memcpy(&bits, &x[k], 8);
    if (bits) {
      idx.push_back(k);
      val.push_back(x[k]);
    }
  }
  Obj o;
  o.vi("i", idx).vf("v", val);
  return o;
}

// Fixed right-hand sides: a dense random vector and a sparse one (every 5th
// position from 1, random values), drawn from HighsRandom(seed).
struct Rhs {
  std::vector<double> dense;
  std::vector<lu_int> sparse_i;
  std::vector<double> sparse_v;
};

Rhs makeRhs(lu_int m, HighsUInt seed) {
  HighsRandom rng(seed);
  Rhs r;
  for (lu_int i = 0; i < m; i++) r.dense.push_back(2 * rng.fraction() - 1);
  for (lu_int i = 1; i < m; i += 5) {
    r.sparse_i.push_back(i);
    r.sparse_v.push_back(2 * rng.fraction() - 1);
  }
  return r;
}

// Dense 'N','T' then sparse 'N','T' solves of the fixed RHS; each solve
// records its status and result, the last one also the scalar block (the
// solve counters accumulate, so a diverging solve still shows there).
std::string solveBatch(struct basiclu_object* obj, lu_int m, const Rhs& rhs) {
  std::string arr = "[";
  bool first = true;
  const char trans[2] = {'N', 'T'};
  for (char t : trans) {
    std::vector<double> lhs(m, 0.0);
    lu_int status = basiclu_obj_solve_dense(obj, rhs.dense.data(), lhs.data(), t);
    Obj o;
    o.s("op", "dense").s("trans", std::string(1, t)).i("status", status)
        .vf("lhs", lhs);
    if (!first) arr += ',';
    first = false;
    arr += o.str();
  }
  for (char t : trans) {
    lu_int status = basiclu_obj_solve_sparse(obj, (lu_int)rhs.sparse_i.size(),
                                             rhs.sparse_i.data(), rhs.sparse_v.data(), t);
    Obj o;
    o.s("op", "sparse").s("trans", std::string(1, t)).i("status", status)
        .i("nzlhs", obj->nzlhs).vi("ilhs", obj->ilhs, obj->nzlhs > 0 ? obj->nzlhs : 0)
        .o("lhs", nzObj(obj->lhs, m));
    if (t == 'T') o.o("xhead", xheadObj(obj->xstore));
    arr += ',';
    arr += o.str();
  }
  arr += ']';
  return arr;
}

// Entering column of variable q of [A I].
void enteringColumn(const HighsSparseMatrix& a, HighsInt q, std::vector<lu_int>& col_i,
                    std::vector<double>& col_v) {
  col_i.clear();
  col_v.clear();
  if (q < a.num_col_) {
    for (HighsInt e = a.start_[q]; e < a.start_[q + 1]; e++) {
      col_i.push_back(a.index_[e]);
      col_v.push_back(a.value_[e]);
    }
  } else {
    col_i.push_back(q - a.num_col_);
    col_v.push_back(1.0);
  }
}

template <typename T>
T* dupArray(const T* src, size_t n) {
  T* dst = (T*)malloc((n ? n : 1) * sizeof(T));
  if (n) memcpy(dst, src, n * sizeof(T));
  return dst;
}

// Deep copy of a BASICLU object (scratch for candidate screening).
void objCopy(struct basiclu_object* dst, const struct basiclu_object* src, lu_int m) {
  const double* xs = src->xstore;
  const size_t lmem = (size_t)xs[BASICLU_MEMORYL], umem = (size_t)xs[BASICLU_MEMORYU],
               wmem = (size_t)xs[BASICLU_MEMORYW];
  dst->istore = dupArray(src->istore, BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * m);
  dst->xstore = dupArray(src->xstore, BASICLU_SIZE_XSTORE_1 + BASICLU_SIZE_XSTORE_M * m);
  dst->Li = dupArray(src->Li, lmem);
  dst->Lx = dupArray(src->Lx, lmem);
  dst->Ui = dupArray(src->Ui, umem);
  dst->Ux = dupArray(src->Ux, umem);
  dst->Wi = dupArray(src->Wi, wmem);
  dst->Wx = dupArray(src->Wx, wmem);
  dst->lhs = dupArray(src->lhs, (size_t)m);
  dst->ilhs = dupArray(src->ilhs, (size_t)m);
  dst->nzlhs = src->nzlhs;
  dst->realloc_factor = src->realloc_factor;
}

Obj inputObj(const BasisMatrix& b) {
  Obj o;
  o.vi("bbegin", b.begin).vi("bend", b.end).vi("bi", b.index).vf("bx", b.value);
  return o;
}

void factorizeCase(Out& out, const std::string& name, const std::string& kind,
                   const HighsSparseMatrix& a, std::vector<HighsInt> basic_index,
                   HighsInt num_updates, double abstol = -1, HighsInt solve_every = 1) {
  const lu_int m = a.num_row_;
  const HighsInt n_tot = a.num_col_ + a.num_row_;
  const Rhs rhs = makeRhs(m, 11);
  BasisMatrix b = basisMatrix(a, basic_index);
  struct basiclu_object obj;
  objInitZeroed(&obj, m);
  // abstol >= 0: IPX's strict mode (BasicLu::_Factorize with
  // strict_abs_pivottol sets BASICLU_ABS_PIVOT_TOLERANCE = kLuDependencyTol).
  if (abstol >= 0) obj.xstore[BASICLU_ABS_PIVOT_TOLERANCE] = abstol;
  lu_int status = basiclu_obj_factorize(&obj, b.begin.data(), b.end.data(),
                                        b.index.data(), b.value.data());
  out.put(Obj("factorize").s("lp", name).s("basis", kind).i("m", m)
              .vi("basic_index", basic_index).o("b", inputObj(b))
              .f("abstol", obj.xstore[BASICLU_ABS_PIVOT_TOLERANCE]).vf("rhs_dense", rhs.dense).vi("rhs_sparse_i", rhs.sparse_i)
              .vf("rhs_sparse_v", rhs.sparse_v).i("status", status)
              .i("lmem", (lu_int)obj.xstore[BASICLU_MEMORYL])
              .i("umem", (lu_int)obj.xstore[BASICLU_MEMORYU])
              .i("wmem", (lu_int)obj.xstore[BASICLU_MEMORYW])
              .o("xhead", xheadObj(obj.xstore)).o("factors", factorsObj(&obj, m))
              .raw("solves", solveBatch(&obj, m, rhs)));
  if (status != BASICLU_OK) {
    basiclu_obj_free(&obj);
    return;
  }

  // Forrest-Tomlin sequence. Entering q: random nonbasic variable; leaving
  // position p = first maximum of |lhs| over ilhs[0..nzlhs) (pattern order)
  // of its tableau column, rejected below 1e-7. Candidates are screened on a
  // scratch copy of the object, so the recorded object sees exactly one
  // solve_for_update 'N' (q), one 'T' (p) and update(xtbl = lhs[p]) per step.
  std::vector<int8_t> is_basic(n_tot, 0);
  for (HighsInt v : basic_index) is_basic[v] = 1;
  HighsRandom rng(2000 + m);
  for (HighsInt k = 0; k < num_updates; k++) {
    HighsInt q = -1;
    lu_int p = -1;
    std::vector<lu_int> col_i;
    std::vector<double> col_v;
    for (HighsInt tries = 0; tries < 200 && q < 0; tries++) {
      HighsInt cand = rng.integer(n_tot);
      if (is_basic[cand]) continue;
      enteringColumn(a, cand, col_i, col_v);
      struct basiclu_object scratch;
      objCopy(&scratch, &obj, m);
      lu_int st = basiclu_obj_solve_for_update(&scratch, (lu_int)col_i.size(), col_i.data(),
                                               col_v.data(), 'N', 1);
      if (st == BASICLU_OK) {
        double best = 0;
        lu_int best_p = -1;
        for (lu_int e = 0; e < scratch.nzlhs; e++) {
          lu_int r = scratch.ilhs[e];
          double v = std::fabs(scratch.lhs[r]);
          if (v > best) {
            best = v;
            best_p = r;
          }
        }
        if (best_p >= 0 && best >= 1e-7) {
          q = cand;
          p = best_p;
        }
      }
      basiclu_obj_free(&scratch);
    }
    if (q < 0) break;
    enteringColumn(a, q, col_i, col_v);
    lu_int ftran_status = basiclu_obj_solve_for_update(&obj, (lu_int)col_i.size(),
                                                       col_i.data(), col_v.data(), 'N', 1);
    const double xtbl = obj.lhs[p];
    Obj o("update_step");
    o.s("lp", name).s("basis", kind).i("k", k).i("variable_in", q).i("pos_out", p)
        .i("variable_out", basic_index[p]).vi("col_i", col_i).vf("col_v", col_v)
        .i("ftran_status", ftran_status).i("ftran_nzlhs", obj.nzlhs)
        .vi("ftran_ilhs", obj.ilhs, obj.nzlhs > 0 ? obj.nzlhs : 0).o("ftran_lhs", nzObj(obj.lhs, m));
    lu_int btran_status = ftran_status == BASICLU_OK
                              ? basiclu_obj_solve_for_update(&obj, 0, &p, NULL, 'T', 1)
                              : ftran_status;
    o.i("btran_status", btran_status).i("btran_nzlhs", obj.nzlhs)
        .vi("btran_ilhs", obj.ilhs, obj.nzlhs > 0 ? obj.nzlhs : 0).o("btran_lhs", nzObj(obj.lhs, m));
    lu_int update_status = btran_status == BASICLU_OK ? basiclu_obj_update(&obj, xtbl) : btran_status;
    o.f("xtbl", xtbl).i("update_status", update_status)
        .i("lmem", (lu_int)obj.xstore[BASICLU_MEMORYL])
        .i("umem", (lu_int)obj.xstore[BASICLU_MEMORYU])
        .i("wmem", (lu_int)obj.xstore[BASICLU_MEMORYW])
        .o("xhead", xheadObj(obj.xstore));
    // Long sequences record the solve batch only every solve_every steps
    // (and on the last one) to keep the fixture small.
    if (update_status == BASICLU_OK && ((k + 1) % solve_every == 0 || k + 1 == num_updates))
      o.raw("solves", solveBatch(&obj, m, rhs));
    out.put(o);
    if (update_status != BASICLU_OK) break;
    is_basic[basic_index[p]] = 0;
    is_basic[q] = 1;
    basic_index[p] = q;
  }

  // Refactorize the updated basis on the used store.
  BasisMatrix b2 = basisMatrix(a, basic_index);
  status = basiclu_obj_factorize(&obj, b2.begin.data(), b2.end.data(), b2.index.data(),
                                 b2.value.data());
  out.put(Obj("refactorize").s("lp", name).s("basis", kind).vi("basic_index", basic_index)
              .o("b", inputObj(b2)).i("status", status)
              .i("lmem", (lu_int)obj.xstore[BASICLU_MEMORYL])
              .i("umem", (lu_int)obj.xstore[BASICLU_MEMORYU])
              .i("wmem", (lu_int)obj.xstore[BASICLU_MEMORYW])
              .o("xhead", xheadObj(obj.xstore)).o("factors", factorsObj(&obj, m))
              .raw("solves", solveBatch(&obj, m, rhs)));
  basiclu_obj_free(&obj);
}

}  // namespace

int runBasiclu(const std::string& mps, const std::string& name, const std::string& path) {
  EkkSnapshot opt = oracleSnapshot(mps, -1);
  if (!opt.valid) {
    fprintf(stderr, "highs-oracle basiclu: no HEkk snapshot for %s\n", mps.c_str());
    return 1;
  }
  const HighsSparseMatrix& a = opt.scaled_lp.a_matrix_;
  const HighsInt num_col = a.num_col_, num_row = a.num_row_;
  EkkSnapshot mid = oracleSnapshot(mps, opt.iteration_count / 2);

  Out out(path);
  out.put(Obj("basiclu_lp").s("name", name).i("num_col", num_col).i("num_row", num_row)
              .i("optimal_iteration_count", opt.iteration_count)
              .i("mid_iteration_count", mid.iteration_count));

  std::vector<HighsInt> slack(num_row);
  for (HighsInt i = 0; i < num_row; i++) slack[i] = num_col + i;
  // Rank deficient: the optimal basis with every third basic variable
  // replaced by the first structural basic variable (factor.cpp's "dup").
  std::vector<HighsInt> dup = opt.basic_index;
  HighsInt first_struct = -1;
  for (HighsInt v : dup)
    if (v < num_col) {
      first_struct = v;
      break;
    }
  if (first_struct < 0) first_struct = 0;
  for (HighsInt i = 1; i < num_row; i += 3) dup[i] = first_struct;

  const HighsInt num_updates = 10;
  factorizeCase(out, name, "slack", a, slack, num_updates);
  factorizeCase(out, name, "mid", a, mid.basic_index, num_updates);
  factorizeCase(out, name, "optimal", a, opt.basic_index, num_updates);
  factorizeCase(out, name, "rank_deficient_dup", a, dup, 0);
  // kLuDependencyTol (ipx_internal.h) = 1e-3.
  factorizeCase(out, name, "optimal_strict", a, opt.basic_index, 4, 1e-3);
  factorizeCase(out, name, "rank_deficient_dup_strict", a, dup, 0, 1e-3);
  // Long Forrest-Tomlin run (up to m-1 updates): reaches the permute reset
  // (pivotlen + nreach > 2m) and the L/U compression in lu_update.
  factorizeCase(out, name, "slack_long", a, slack, a.num_row_ - 1, -1, 8);
  return 0;
}

// Synthetic matrices for the pivot routines the small Netlib bases do not
// reach: a dense m x 2m matrix (pivot columns with more than MAXROW_SMALL =
// 64 off-diagonals take lu_pivot_any), a random sparse one, and a dense
// 40 x 80 one whose basis columns m-4..m-1 are columns 0..3 plus 1e-6 noise
// (neardep40: after elimination their maximum drops below the absolute pivot
// tolerance, so lu_pivot's lu_remove_col and the rank-deficient Markowitz
// exits run). Basis = the first m columns; entering columns are drawn from
// all 3m variables of [A I] as for the LPs.
int runBasicluSynth(const std::string& name, const std::string& path) {
  HighsInt m;
  double density;
  if (name == "dense70") {
    m = 70;
    density = 1.0;
  } else if (name == "rand120") {
    m = 120;
    density = 0.06;
  } else if (name == "neardep40") {
    m = 40;
    density = 1.0;
  } else {
    fprintf(stderr, "highs-oracle basiclu-synth: unknown matrix %s\n", name.c_str());
    return 64;
  }
  HighsSparseMatrix a;
  a.format_ = MatrixFormat::kColwise;
  a.num_row_ = m;
  a.num_col_ = 2 * m;
  HighsRandom rng(31 + m);
  a.start_.assign(1, 0);  // HighsSparseMatrix() already holds start_ = {0}
  a.index_.clear();
  a.value_.clear();
  for (HighsInt j = 0; j < a.num_col_; j++) {
    for (HighsInt i = 0; i < m; i++) {
      if (density < 1.0 && rng.fraction() >= density && i != j % m) continue;
      a.index_.push_back(i);
      a.value_.push_back(2 * rng.fraction() - 1);
    }
    a.start_.push_back((HighsInt)a.index_.size());
  }
  if (name == "neardep40") {
    for (HighsInt k = 0; k < 4; k++) {
      const HighsInt j = m - 4 + k;
      for (HighsInt i = 0; i < m; i++)
        a.value_[a.start_[j] + i] = a.value_[a.start_[k] + i] + 1e-6 * (2 * rng.fraction() - 1);
    }
  }
  std::vector<HighsInt> first(m);
  for (HighsInt i = 0; i < m; i++) first[i] = i;
  Out out(path);
  out.put(Obj("basiclu_lp").s("name", name).i("num_col", a.num_col_).i("num_row", m)
              .f("density", density));
  if (name == "neardep40") {
    factorizeCase(out, name, "first_m", a, first, 3);
    factorizeCase(out, name, "first_m_strict", a, first, 3, 1e-3);
    return 0;
  }
  factorizeCase(out, name, "first_m", a, first, 6);
  return 0;
}

