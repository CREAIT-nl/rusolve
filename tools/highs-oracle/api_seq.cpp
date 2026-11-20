// api_seq: the C++ half of the rusolve-ekk model-modification parity
// harness. It drives the
// real HiGHS public API through the call sequence EkkEngine uses and prints
// one line per stage; `rusolve-cli ekk-seq <file>` prints the same lines from
// the Rust port, and the two outputs must be identical.
//
// Sequence (all on one Highs object, presolve off after stage 0, as the MIP
// LP relaxation runs it):
//   0 readModel + clearIntegrality + run          (cold, default presolve)
//   1 changeColsBounds(set) on one column + run   (hot start: kNewBounds)
//   2 addRows(one row) + run                      (HEkk cleared, basis kept)
//   3 deleteRows(mask) of that row + run          (HEkk cleared, basis compacted)
// Stage 1 column: the first j with the largest col_value (ties -> lowest j);
// new upper = max(col_lower[j], 0.5 * col_value[j]) when col_value[j] > 0,
// otherwise new lower = min(col_upper[j], col_value[j] + 1).
// Stage 2 row: coefficient 1.0 on every column with col_value != 0 after
// stage 1 (ascending j), lower = -inf, upper = 0.9 * sum of those values
// when the sum is > 0, else sum - 1. No such column: the row is skipped
// (the line prints "skip").
// Output line: "stage <k> status=<model status> iters=<n> obj=<%.17g>".
//
// Build (against the -ffp-contract=off oracle build of HiGHS):
//   clang++ -std=c++11 -O2 -ffp-contract=off tools/highs-oracle/api_seq.cpp \
//     -I target/highs-oracle/highs-src/highs \
//     -I target/highs-oracle/highs-src/extern \
//     -I target/highs-oracle/build/highs-build \
//     target/highs-oracle/build/lib/libhighs.a \
//     target/highs-oracle/build/lib/libhighs_extras.a \
//     -o target/highs-oracle/api_seq

// Or, portably (gen_fixtures.sh builds it into target/highs-oracle/build/):
//   cmake --build target/highs-oracle/build --target mip_trace api_seq

#include <cmath>
#include <cstdio>
#include <vector>

#include "Highs.h"

static void report(Highs& h, int stage) {
  std::printf("stage %d status=%s iters=%d obj=%.17g\n", stage,
              h.modelStatusToString(h.getModelStatus()).c_str(),
              (int)h.getInfo().simplex_iteration_count,
              h.getInfo().objective_function_value);
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: api_seq <model.mps>\n");
    return 2;
  }
  Highs h;
  h.setOptionValue("output_flag", false);
  if (h.readModel(argv[1]) == HighsStatus::kError) return 1;
  h.clearIntegrality();
  h.run();
  report(h, 0);
  h.setOptionValue("presolve", "off");

  const HighsLp& lp = h.getLp();
  const HighsInt n = lp.num_col_;
  if (n == 0) return 0;
  std::vector<double> x = h.getSolution().col_value;
  HighsInt j = 0;
  for (HighsInt k = 1; k < n; k++)
    if (x[k] > x[j]) j = k;
  double lo = lp.col_lower_[j], up = lp.col_upper_[j];
  if (x[j] > 0)
    up = std::max(lo, 0.5 * x[j]);
  else
    lo = std::min(up, x[j] + 1);
  h.changeColsBounds(1, &j, &lo, &up);
  h.run();
  report(h, 1);

  x = h.getSolution().col_value;
  std::vector<HighsInt> index;
  std::vector<double> value;
  double sum = 0;
  for (HighsInt k = 0; k < n; k++)
    if (x[k] != 0) {
      index.push_back(k);
      value.push_back(1.0);
      sum += x[k];
    }
  if (index.empty()) {
    std::printf("stage 2 skip\n");
    return 0;
  }
  const double row_lower = -kHighsInf;
  const double row_upper = sum > 0 ? 0.9 * sum : sum - 1;
  const HighsInt start = 0;
  h.addRows(1, &row_lower, &row_upper, (HighsInt)index.size(), &start,
            index.data(), value.data());
  h.run();
  report(h, 2);

  std::vector<HighsInt> mask(h.getNumRow(), 0);
  mask.back() = 1;
  h.deleteRows(mask.data());
  h.run();
  report(h, 3);
  return 0;
}
