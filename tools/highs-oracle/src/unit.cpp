// Unit-level fixtures: HighsSparseMatrix PRICE / collectAj / partitioned
// row-wise matrix maintenance, and HVector primitives, on a real LP matrix.
// Every record embeds its inputs, so a Rust test can replay it alone.
#include <cmath>

#include "lp_data/HConst.h"
#include "oracle.h"
#include "util/HighsRandom.h"

namespace {

HVector makeVec(HighsInt size) {
  HVector v;
  v.setup(size);
  v.clear();
  return v;
}

// Deterministic RHS-like vectors of dimension n: unit, sparse (3 entries),
// medium (~n/4 entries), dense (every entry). Values in (-1,1) from
// HighsRandom(seed); a few exact cancellations are not attempted here.
std::vector<std::pair<std::string, HVector>> inputVectors(HighsInt n,
                                                          HighsUInt seed) {
  HighsRandom rng(seed);
  std::vector<std::pair<std::string, HVector>> out;
  {
    HVector v = makeVec(n);
    HighsInt r = n / 2;
    v.index[v.count++] = r;
    v.array[r] = 1.0;
    out.emplace_back("unit", v);
  }
  {
    HVector v = makeVec(n);
    HighsInt picks[3] = {n - 1, 0, n / 3};
    for (HighsInt r : picks) {
      if (v.array[r] != 0) continue;
      v.index[v.count++] = r;
      v.array[r] = 2 * rng.fraction() - 1;
    }
    out.emplace_back("sparse", v);
  }
  {
    HVector v = makeVec(n);
    for (HighsInt r = n - 1; r >= 0; r -= 4) {
      v.index[v.count++] = r;
      v.array[r] = 2 * rng.fraction() - 1;
    }
    out.emplace_back("medium", v);
  }
  {
    HVector v = makeVec(n);
    for (HighsInt r = 0; r < n; r++) {
      v.index[v.count++] = r;
      v.array[r] = 2 * rng.fraction() - 1;
    }
    out.emplace_back("dense", v);
  }
  return out;
}

}  // namespace

int runUnit(const std::string& mps, const std::string& name,
            const std::string& path) {
  HighsLp lp = oracleLoadLp(mps);
  const HighsSparseMatrix& a = lp.a_matrix_;
  const HighsInt num_col = lp.num_col_, num_row = lp.num_row_;
  Out out(path);
  out.put(Obj("lp").s("name", name).i("num_col", num_col).i("num_row", num_row)
              .mat("a", a));

  // createRowwise
  HighsSparseMatrix ar;
  ar.createRowwise(a);
  out.put(Obj("createRowwise").mat("out", ar));

  // createRowwisePartitioned: all-in-partition (NULL), and a mixed partition.
  HighsSparseMatrix arp_all;
  arp_all.createRowwisePartitioned(a, NULL);
  out.put(Obj("createRowwisePartitioned").b("null_partition", true)
              .mat("out", arp_all));
  std::vector<int8_t> part(num_col + num_row, 1);
  for (HighsInt j = 0; j < num_col; j++) part[j] = (j % 3 == 1) ? 0 : 1;
  for (HighsInt i = 0; i < num_row; i++) part[num_col + i] = (i % 3 == 1) ? 1 : 0;
  HighsSparseMatrix arp;
  arp.createRowwisePartitioned(a, part.data());
  out.put(Obj("createRowwisePartitioned").b("null_partition", false)
              .vi8("in_partition", part.data(), num_col).mat("out", arp));

  // PRICE on row_ep-like vectors.
  for (auto& named : inputVectors(num_row, 20261030)) {
    const HVector& col = named.second;
    for (int quad = 0; quad <= 1; quad++) {
      HVector res = makeVec(num_col);
      a.priceByColumn(quad, res, col);
      out.put(Obj("priceByColumn").s("in_kind", named.first).i("quad", quad)
                  .hv("in", col).hv("out", res));
      HVector res_r = makeVec(num_col);
      ar.priceByRow(quad, res_r, col);
      out.put(Obj("priceByRow").s("matrix", "rowwise").s("in_kind", named.first)
                  .i("quad", quad).hv("in", col).hv("out", res_r));
      HVector res_p = makeVec(num_col);
      arp.priceByRow(quad, res_p, col);
      out.put(Obj("priceByRow").s("matrix", "partitioned")
                  .s("in_kind", named.first).i("quad", quad).hv("in", col)
                  .hv("out", res_p));
      // (expected_density, switch_density, from_index): the dual simplex
      // passes switch_density = kHyperPriceDensity and from_index = 0; the
      // other combinations exercise the mid-way switch and a non-zero start.
      struct Combo {
        double ed, sd;
        HighsInt from;
      };
      const Combo combos[4] = {{0.01, kHyperPriceDensity, 0},
                               {0.5, kHyperPriceDensity, 0},
                               {0.01, kHyperPriceDensity, 1},
                               {0.01, 0.02, 0}};
      for (const Combo& c : combos) {
        const double ed = c.ed, sd = c.sd;
        const HighsInt from = c.from;
        {
          {
            if (from >= col.count) continue;
            for (int which = 0; which < 2; which++) {
              const HighsSparseMatrix& m = which ? arp : ar;
              HVector res_s = makeVec(num_col);
              m.priceByRowWithSwitch(quad, res_s, col, ed, from, sd);
              out.put(Obj("priceByRowWithSwitch")
                          .s("matrix", which ? "partitioned" : "rowwise")
                          .s("in_kind", named.first).i("quad", quad)
                          .f("expected_density", ed).i("from_index", from)
                          .f("switch_density", sd).hv("in", col)
                          .hv("out", res_s));
            }
          }
        }
      }
    }
    // computeDot on every structural and three logical columns
    std::vector<double> dots;
    std::vector<HighsInt> use;
    for (HighsInt j = 0; j < num_col; j++) use.push_back(j);
    use.push_back(num_col);
    use.push_back(num_col + num_row / 2);
    use.push_back(num_col + num_row - 1);
    for (HighsInt j : use) dots.push_back(a.computeDot(col, j));
    out.put(Obj("computeDot").s("in_kind", named.first).hv("in", col)
                .vi("use_col", use).vf("out", dots));
  }

  // collectAj: accumulate a sequence of (column, multiplier) into one vector,
  // including an exact cancellation (-> kHighsZero entries) and logicals.
  {
    HighsRandom rng(7);
    HVector v = makeVec(num_row);
    std::vector<HighsInt> seq_col;
    std::vector<double> seq_mul;
    for (HighsInt k = 0; k < 6; k++) {
      HighsInt j = rng.integer(num_col + num_row);
      double mul = 4 * rng.fraction() - 2;
      seq_col.push_back(j);
      seq_mul.push_back(mul);
    }
    // exact cancellation of the first column
    seq_col.push_back(seq_col[0]);
    seq_mul.push_back(-seq_mul[0]);
    // a logical twice
    seq_col.push_back(num_col);
    seq_mul.push_back(1.5);
    seq_col.push_back(num_col);
    seq_mul.push_back(-1.5);
    Obj steps;
    std::string arr = "[";
    for (size_t k = 0; k < seq_col.size(); k++) {
      a.collectAj(v, seq_col[k], seq_mul[k]);
      Obj st;
      st.i("col", seq_col[k]).f("multiplier", seq_mul[k]).hv("after", v);
      if (k) arr += ',';
      arr += st.str();
    }
    arr += ']';
    out.put(Obj("collectAj").i("size", num_row).raw("steps", arr));
  }

  // HighsSparseMatrix::update on the partitioned matrix: in_partition plays
  // the role of nonbasicFlag (1 = nonbasic, in [start, p_end)).
  {
    HighsSparseMatrix m = arp;
    std::vector<int8_t> flag = part;
    HighsRandom rng(11);
    std::string arr = "[";
    for (HighsInt k = 0; k < 8; k++) {
      HighsInt var_in = -1, var_out = -1;
      for (HighsInt tries = 0; tries < 1000 && var_in < 0; tries++) {
        HighsInt j = rng.integer(num_col + num_row);
        if (flag[j]) var_in = j;
      }
      for (HighsInt tries = 0; tries < 1000 && var_out < 0; tries++) {
        HighsInt j = rng.integer(num_col + num_row);
        if (!flag[j]) var_out = j;
      }
      if (var_in < 0 || var_out < 0) break;
      m.update(var_in, var_out, a);
      flag[var_in] = 0;
      flag[var_out] = 1;
      Obj st;
      st.i("var_in", var_in).i("var_out", var_out).mat("after", m);
      if (k) arr += ',';
      arr += st.str();
    }
    arr += ']';
    out.put(Obj("updatePartitioned").vi8("in_partition", part.data(), num_col)
                .mat("before", arp).raw("steps", arr));
  }

  // HVector primitives: tight, pack, reIndex, norm2, saxpy, copy.
  {
    HighsRandom rng(3);
    HVector v = makeVec(num_row);
    for (HighsInt r = 0; r < num_row; r += 2) {
      v.index[v.count++] = r;
      double x = 2 * rng.fraction() - 1;
      if (r % 6 == 0) x *= 1e-14;  // straddles kHighsTiny
      if (r % 10 == 0) x = -0.0;
      v.array[r] = x;
    }
    HVector t = v;
    t.tight();
    out.put(Obj("hvector_tight").hv("in", v).hv("out", t));
    HVector tn = v;
    tn.count = -1;
    tn.tight();
    out.put(Obj("hvector_tight").hv("in", [&] { HVector c = v; c.count = -1; return c; }())
                .hv("out", tn));
    HVector p = t;
    p.packFlag = true;
    p.pack();
    out.put(Obj("hvector_pack").hv("in", t).hv("out", p, true));
    out.put(Obj("hvector_norm2").hv("in", t).f("out", t.norm2()));
    HVector ri = t;
    ri.reIndex();
    out.put(Obj("hvector_reIndex").hv("in", t).hv("out", ri));
    HVector ri2 = t;
    ri2.count = -1;
    ri2.reIndex();
    out.put(Obj("hvector_reIndex").hv("in", [&] { HVector c = t; c.count = -1; return c; }())
                .hv("out", ri2));
    for (auto& named : inputVectors(num_row, 99)) {
      HVector s = t;
      const double mult = -0.75;
      s.saxpy(mult, &named.second);
      out.put(Obj("hvector_saxpy").f("pivotX", mult).hv("in", t)
                  .hv("pivot", named.second).hv("out", s));
    }
  }
  return 0;
}

int runRandom(const std::string& path) {
  Out out(path);
  const HighsUInt seeds[4] = {0, 1, 42, 2147483647u};
  for (HighsUInt seed : seeds) {
    HighsRandom rng(seed);
    std::vector<HighsInt> ints;
    for (int k = 0; k < 16; k++) ints.push_back(rng.integer());
    const HighsInt sups[8] = {1, 2, 3, 7, 100, 1000, 1 << 20, 2147483647};
    std::vector<HighsInt> sup_v, bounded;
    for (int k = 0; k < 64; k++) {
      HighsInt sup = sups[k % 8];
      sup_v.push_back(sup);
      bounded.push_back(rng.integer(sup));
    }
    std::vector<double> frac, cfrac, reals;
    for (int k = 0; k < 16; k++) frac.push_back(rng.fraction());
    for (int k = 0; k < 16; k++) cfrac.push_back(rng.closedFraction());
    for (int k = 0; k < 16; k++) reals.push_back(rng.real(-3.0, 5.0));
    std::vector<HighsInt> bits;
    for (int k = 0; k < 16; k++) bits.push_back(rng.bit());
    std::vector<HighsInt> perm;
    for (HighsInt k = 0; k < 20; k++) perm.push_back(k);
    rng.shuffle(perm.data(), (HighsInt)perm.size());
    // The order of draws above is part of the fixture: integer() x16,
    // integer(sup) x64, fraction() x16, closedFraction() x16,
    // real(-3,5) x16, bit() x16, shuffle(0..20).
    out.put(Obj("random").i("seed", seed).vi("integer", ints)
                .vi("sup", sup_v).vi("integer_sup", bounded).vf("fraction", frac)
                .vf("closedFraction", cfrac).vf("real_m3_5", reals)
                .vi("bit", bits).vi("shuffle20", perm));
  }
  return 0;
}
