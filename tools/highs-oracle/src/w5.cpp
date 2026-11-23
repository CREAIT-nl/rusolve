// W5 oracle commands: the existing simplex and HiPO traces at threads = N with
// parallel = on. Records: FORMAT.md "W5 records" (the W1 simplex / W4 HiPO
// record formats, with "threads" and "parallel" in the header).
//
//   highs-oracle w5-simplex <file.mps> <name> <2|3> <threads> <out.jsonl>
//   highs-oracle w5-hipo    <file.mps> <name> <tree|node|both>
//                           <choose|augmented|normaleq> <threads> <out.jsonl>
//
// w5-simplex is runSimplex (simplex.cpp) with simplex_strategy 2 (SIP) or 3
// (PAMI), threads=N and parallel=on; simplex_max_concurrency stays at its
// default. The simplex hooks fire on the main thread, so the trace does not
// depend on the scheduling.
//
// w5-hipo is the W4 run (w4.cpp, HIPO=ON build on refblas) with threads=N,
// parallel=on, hipo_parallel_type and hipo_system. hipo_system=choose is
// accepted only at threads=1: chooseNla runs the AS analyse and the NE
// structure as two tasks, which race at threads > 1. The hipo_ord records are
// written in ordering order from the calling thread (highs-oracle-w5-hooks.patch).
#include <cstdio>
#include <cstdlib>
#include <string>

#include "oracle.h"

namespace {

// A positive thread count, or 0 on a malformed argument.
int w5Threads(const std::string& s) {
  char* end = nullptr;
  const long v = std::strtol(s.c_str(), &end, 10);
  if (s.empty() || *end != '\0' || v < 1 || v > 64) return 0;
  return (int)v;
}

}  // namespace

int runW5Simplex(const std::string& mps, const std::string& name, const std::string& strategy,
                 const std::string& threads, const std::string& path) {
  if (strategy != "2" && strategy != "3") {
    fprintf(stderr, "highs-oracle: w5-simplex: strategy must be 2 (SIP) or 3 (PAMI)\n");
    return 64;
  }
  const int n = w5Threads(threads);
  if (n == 0) {
    fprintf(stderr, "highs-oracle: w5-simplex: bad thread count %s\n", threads.c_str());
    return 64;
  }
  return runSimplexThreads(mps, name, atoi(strategy.c_str()), n, path);
}

int runW5Hipo(const std::string& mps, const std::string& name, const std::string& type,
              const std::string& system, const std::string& threads, const std::string& path) {
  if (type != "tree" && type != "node" && type != "both") {
    fprintf(stderr, "highs-oracle: w5-hipo: parallel type must be tree, node or both\n");
    return 64;
  }
  if (system != "choose" && system != "augmented" && system != "normaleq") {
    fprintf(stderr, "highs-oracle: w5-hipo: system must be choose, augmented or normaleq\n");
    return 64;
  }
  const int n = w5Threads(threads);
  if (n == 0) {
    fprintf(stderr, "highs-oracle: w5-hipo: bad thread count %s\n", threads.c_str());
    return 64;
  }
  if (system == "choose" && n != 1) {
    fprintf(stderr,
            "highs-oracle: w5-hipo: hipo_system=choose races at threads > 1 "
            "(its two analyse tasks race); use augmented or normaleq\n");
    return 64;
  }
  return runW4WithOverrides(mps, name, path,
                            "solver=hipo,threads=" + std::to_string(n) +
                                ",parallel=on,hipo_parallel_type=" + type +
                                ",hipo_system=" + system,
                            n);
}
