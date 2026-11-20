// Clique-row normalisation oracle (1.16.0 form of HiGHS #3286 + #3318): one call of
// HPresolve::normaliseCliqueRows on a model, the presolve state after it.
//
//   highs-oracle g1-clique-normalise <model> <name> default <out.jsonl>
//
// The set-up is the one of the 1.16.0 TestPresolveRules normalisation cases
// (check/TestPresolveRules.cpp, "test-normalise-unequal-coeff"): a
// HighsMipSolver on the model with mipdata_ init and
// setupDomainPropagation, then HPresolve::setInput(mipsolver, -1),
// okSetupPresolveDataStructures and normaliseCliqueRows. Records:
// g1_header (as every G1 command) and clqn_result (FORMAT.md "G1 records").
//
// normaliseCliqueRows does not exist before 1.16 (#3286 is after 1.15.1):
// on the 1.15.1 oracle the command fails with 64 and gen_fixtures.sh only
// runs it with ORACLE_HIGHS=1.16.0.
//
// The presolve state is read with `#define private public` and, since most
// HPresolve members precede any access specifier, `#define class struct`
// (this translation unit only, around HPresolve.h alone, after every standard
// and HiGHS header it includes; the header's only other uses of the keyword
// are a forward declaration and `enum class`, both valid as `struct`).
#include <cmath>
#include <cstdio>
#include <list>
#include <memory>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "HConfig.h"
#include "Highs.h"
#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsCallback.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "oracle.h"
#include "parallel/HighsParallel.h"
#include "presolve/HPresolveAnalysis.h"
#if !(HIGHS_VERSION_MAJOR == 1 && HIGHS_VERSION_MINOR < 16)
#include "presolve/HPresolveCliqueTable.h"
#endif
#include "presolve/HighsPostsolveStack.h"
#include "util/HighsCDouble.h"
#include "util/HighsHash.h"
#include "util/HighsHashTree.h"
#include "util/HighsLinearSumBounds.h"
#include "util/HighsMatrixSlice.h"

#define private public
#define protected public
#define class struct
#include "presolve/HPresolve.h"
#undef class
#undef private
#undef protected

int runClqn(const std::string& model, const std::string& name, const std::string& scenario,
            const std::string& path) {
#if HIGHS_VERSION_MAJOR == 1 && HIGHS_VERSION_MINOR < 16
  (void)model;
  (void)name;
  (void)scenario;
  (void)path;
  fprintf(stderr, "highs-oracle: g1-clique-normalise needs HiGHS 1.16\n");
  return 64;
#else
  if (scenario != "default") {
    fprintf(stderr, "highs-oracle: unknown scenario %s\n", scenario.c_str());
    return 64;
  }
  const HighsLp lp = oracleG1Read(model);
  Out out(path);
  {
    Obj h("g1_header");
    h.s("name", name).s("highs_version", highsVersion()).s("highs_githash", highsGithash())
        .s("overrides", scenario);
    oracleG1AddLp(h, lp);
    out.put(h);
  }

  highs::parallel::initialize_scheduler(1);
  Highs highs;
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("threads", 1);
  highs.setOptionValue("parallel", "off");
  highs.setOptionValue("random_seed", 0);
  if (highs.passModel(lp) != HighsStatus::kOk) return 2;

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;
  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();
  presolve::HighsPostsolveStack& postsolve_stack = mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  if (!presolve.okSetupPresolveDataStructures()) return 2;
  const presolve::HPresolve::Result result = presolve.normaliseCliqueRows(postsolve_stack);

  // The model after the call: column bounds and deletions, then every row
  // (bounds, deleted flag, nonzeros in getRowVector order, implied activity
  // bounds).
  const HighsLp& m = *presolve.model;
  std::vector<HighsInt> col_deleted, row_deleted, row_start, row_index, num_inf_lower,
      num_inf_upper;
  std::vector<double> row_value, sum_lower, sum_upper;
  for (HighsInt c = 0; c < m.num_col_; c++) col_deleted.push_back(presolve.colDeleted[c] ? 1 : 0);
  for (HighsInt r = 0; r < m.num_row_; r++) {
    row_deleted.push_back(presolve.rowDeleted[r] ? 1 : 0);
    row_start.push_back(static_cast<HighsInt>(row_index.size()));
    if (presolve.rowDeleted[r]) {
      sum_lower.push_back(0.0);
      sum_upper.push_back(0.0);
      num_inf_lower.push_back(0);
      num_inf_upper.push_back(0);
      continue;
    }
    for (const HighsSliceNonzero& nz : presolve.getRowVector(r)) {
      row_index.push_back(nz.index());
      row_value.push_back(nz.value());
    }
    sum_lower.push_back(presolve.impliedRowBounds.getSumLower(r));
    sum_upper.push_back(presolve.impliedRowBounds.getSumUpper(r));
    num_inf_lower.push_back(presolve.impliedRowBounds.getNumInfSumLower(r));
    num_inf_upper.push_back(presolve.impliedRowBounds.getNumInfSumUpper(r));
  }
  row_start.push_back(static_cast<HighsInt>(row_index.size()));

  Obj c("clqn_result");
  c.i("result", static_cast<int>(result))
      .i("num_reductions", static_cast<long long>(postsolve_stack.numReductions()))
      .vf("col_lower", m.col_lower_).vf("col_upper", m.col_upper_).vi("col_deleted", col_deleted)
      .vf("row_lower", m.row_lower_).vf("row_upper", m.row_upper_).vi("row_deleted", row_deleted)
      .vi("row_start", row_start).vi("row_index", row_index).vf("row_value", row_value)
      .vf("sum_lower", sum_lower).vf("sum_upper", sum_upper)
      .vi("num_inf_sum_lower", num_inf_lower).vi("num_inf_sum_upper", num_inf_upper);
  out.put(c);
  HighsTaskExecutor::shutdown(true);
  return 0;
#endif
}
