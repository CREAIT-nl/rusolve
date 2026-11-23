// mip_trace: the C++ half of the rusolve-ekk MIP port parity harness.
// It runs the real HiGHS
// MIP solver serially (threads=1, random_seed=0, default mip_rel_gap) and
// prints a deterministic trace. `rusolve-cli mip-trace <file>` prints the
// same lines from the Rust port; the two outputs must be identical.
//
// Trace lines (all doubles %.17g, so any bit difference shows):
//   log <text>            selected HiGHS log lines (presolve dimensions and
//                         reductions, restarts, the final solving report),
//                         with every "<number>s" timing token replaced by
//                         "<t>s". The time-driven node display lines are
//                         dropped (their cadence depends on wall clock).
//   cb nodes=<n> lpit=<it> db=<dual bound> pb=<primal bound>
//                         one line per kCallbackMipInterrupt call
//                         (HighsMipSolverData::checkLimits: root separation
//                         rounds, then every node of the search) whose
//                         (nodes, db, pb) differs from the previous line.
//   sol nodes=<n> lpit=<it> obj=<objective>
//                         one line per kCallbackMipImprovingSolution.
//   final status=<model status> obj=<objective> bound=<mip_dual_bound>
//         nodes=<mip_node_count> lp_iters=<simplex_iteration_count>
//
// Only the interrupt / improving-solution / logging callbacks are started.
// They read solver state and do not change the search path (the
// kCallbackMipSolution / user-solution / cut-pool callbacks, which do, stay
// off).
//
// Build (against the -ffp-contract=off, serial oracle build of HiGHS):
//   R=target/highs-oracle      # from the checkout root
//   clang++ -std=c++11 -O2 -ffp-contract=off tools/highs-oracle/mip_trace.cpp \
//     -I $R/highs-src/highs -I $R/highs-src/extern -I $R/build/highs-build \
//     $R/build/lib/libhighs.a $R/build/lib/libhighs_extras.a -o $R/mip_trace
// Or, portably (gen_fixtures.sh builds it into target/highs-oracle/build/):
//   cmake --build target/highs-oracle/build --target mip_trace api_seq
// Usage: mip_trace <model.mps> [time_limit_s (default 3600)] [key=value ...]
//
// Each key=value argument (after the model path, and after the optional
// time_limit_s) is applied with Highs::setOptionValue(key, value) on the
// string overload, after the defaults (threads=1, parallel=off, random_seed=0,
// the time limit), so it can override them. An
// unknown option, or a value HiGHS rejects, exits with status 1. Examples:
//   mip_trace X.mps presolve=off mip_max_nodes=1
//   mip_trace X.mps 60 presolve=off
// The option values are not echoed to stdout, so a trace with options can
// be compared line by line with one without.
//
// Analytic-centre and LP-relaxation IPM records (IPX results, which the port
// replays from fixtures) are written to the file named by the environment
// variable HIGHS_ORACLE_AC_DUMP, when it is set. They do not touch stdout
// (see highs-oracle-hooks.patch and FORMAT.md, "ac" / "lp_ipm").

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "Highs.h"

namespace {

struct TraceState {
  bool have_last = false;
  long long last_nodes = 0;
  double last_db = 0;
  double last_pb = 0;
};

bool sameDouble(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

// Replace every "<digits>[.<digits>]s" token by "<t>s".
std::string stripTimes(const std::string& s) {
  std::string out;
  size_t i = 0;
  while (i < s.size()) {
    if (std::isdigit(static_cast<unsigned char>(s[i])) &&
        (i == 0 || !std::isalnum(static_cast<unsigned char>(s[i - 1])))) {
      size_t j = i;
      while (j < s.size() &&
             (std::isdigit(static_cast<unsigned char>(s[j])) || s[j] == '.'))
        ++j;
      if (j < s.size() && s[j] == 's' &&
          (j + 1 == s.size() ||
           !std::isalnum(static_cast<unsigned char>(s[j + 1])))) {
        out += "<t>s";
        i = j + 1;
        continue;
      }
      out.append(s, i, j - i);
      i = j;
      continue;
    }
    out += s[i++];
  }
  return out;
}

bool keepLogLine(const std::string& line) {
  static const char* keys[] = {"Presolv",      "reductions",   "rows,",
                               "Restart",      "Solving report", "Status",
                               "Primal bound", "Dual bound",   "Gap",
                               "Solution status", "Nodes",     "LP iterations",
                               "Objective",    "Coefficient ranges",
                               "Cost", "Bound", "RHS", "Matrix", "MIP has",
                               "symmetr", "Symmetr", "clique", "Clique"};
  if (line.find("Timing") != std::string::npos) return false;
  // the timing block of the solving report
  if (line.find("(Presolve)") != std::string::npos ||
      line.find("(Solve)") != std::string::npos ||
      line.find("(Postsolve)") != std::string::npos)
    return false;
  for (const char* k : keys)
    if (line.find(k) != std::string::npos) return true;
  return false;
}

void callback(int type, const std::string& message,
              const HighsCallbackOutput* out, HighsCallbackInput* in,
              void* user) {
  TraceState& st = *static_cast<TraceState*>(user);
  switch (type) {
    case kCallbackLogging: {
      std::string text = message;
      while (!text.empty() && (text.back() == '\n' || text.back() == '\r'))
        text.pop_back();
      // a log message may hold several lines
      size_t pos = 0;
      while (pos <= text.size()) {
        size_t nl = text.find('\n', pos);
        std::string line =
            text.substr(pos, nl == std::string::npos ? std::string::npos
                                                     : nl - pos);
        if (!line.empty() && keepLogLine(line))
          std::printf("log %s\n", stripTimes(line).c_str());
        if (nl == std::string::npos) break;
        pos = nl + 1;
      }
      break;
    }
    case kCallbackMipInterrupt: {
      long long nodes = static_cast<long long>(out->mip_node_count);
      if (st.have_last && nodes == st.last_nodes &&
          sameDouble(out->mip_dual_bound, st.last_db) &&
          sameDouble(out->mip_primal_bound, st.last_pb))
        break;
      st.have_last = true;
      st.last_nodes = nodes;
      st.last_db = out->mip_dual_bound;
      st.last_pb = out->mip_primal_bound;
      std::printf("cb nodes=%lld lpit=%lld db=%.17g pb=%.17g\n", nodes,
                  static_cast<long long>(out->mip_total_lp_iterations),
                  out->mip_dual_bound, out->mip_primal_bound);
      break;
    }
    case kCallbackMipImprovingSolution:
      std::printf("sol nodes=%lld lpit=%lld obj=%.17g\n",
                  static_cast<long long>(out->mip_node_count),
                  static_cast<long long>(out->mip_total_lp_iterations),
                  out->objective_function_value);
      break;
    default:
      break;
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: mip_trace <model.mps> [time_limit_s] [key=value ...]\n");
    return 2;
  }
  // argv[2] is the time limit when it has no '='; the rest are key=value.
  int first_option = 2;
  double time_limit = 3600.0;
  if (argc > 2 && std::strchr(argv[2], '=') == nullptr) {
    time_limit = std::atof(argv[2]);
    first_option = 3;
  }
  Highs h;
  h.setOptionValue("output_flag", true);
  h.setOptionValue("log_to_console", false);
  h.setOptionValue("threads", 1);
  // HiGHS 1.16 (#3406) runs the MIP in parallel under the default
  // parallel=choose; off keeps the serial search on both oracle versions
  // (1.15.1 already treated choose as off, so its traces are unchanged).
  h.setOptionValue("parallel", "off");
  h.setOptionValue("random_seed", 0);
  h.setOptionValue("time_limit", time_limit);
  for (int i = first_option; i < argc; ++i) {
    const std::string arg = argv[i];
    const size_t eq = arg.find('=');
    if (eq == std::string::npos || eq == 0) {
      std::fprintf(stderr, "mip_trace: expected key=value, got '%s'\n",
                   argv[i]);
      return 2;
    }
    const std::string key = arg.substr(0, eq);
    const std::string value = arg.substr(eq + 1);
    if (h.setOptionValue(key, value) == HighsStatus::kError) {
      std::fprintf(stderr,
                   "mip_trace: unknown option or invalid value: %s=%s\n",
                   key.c_str(), value.c_str());
      return 1;
    }
  }
  if (h.readModel(argv[1]) == HighsStatus::kError) {
    std::fprintf(stderr, "mip_trace: cannot read %s\n", argv[1]);
    return 1;
  }
  TraceState st;
  h.setCallback(callback, &st);
  h.startCallback(kCallbackLogging);
  h.startCallback(kCallbackMipInterrupt);
  h.startCallback(kCallbackMipImprovingSolution);
  h.run();
  const HighsInfo& info = h.getInfo();
  std::printf("final status=%s obj=%.17g bound=%.17g nodes=%lld lp_iters=%lld\n",
              h.modelStatusToString(h.getModelStatus()).c_str(),
              info.objective_function_value, info.mip_dual_bound,
              static_cast<long long>(info.mip_node_count),
              static_cast<long long>(info.simplex_iteration_count));
  return 0;
}
