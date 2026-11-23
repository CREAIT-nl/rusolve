// highs-oracle: golden-fixture generator for the rusolve-ekk faithful port of
// the HiGHS simplex (HiGHS is MIT-licensed, (c) the HiGHS authors; this tool
// links an instrumented copy of it, see ../highs-oracle-hooks.patch).
//
//   highs-oracle unit    <file.mps> <name> <out.jsonl>
//   highs-oracle factor  <file.mps> <name> <out.jsonl>
//   highs-oracle simplex <file.mps> <name> <strategy> <out.jsonl>
//   highs-oracle basiclu <file.mps> <name> <out.jsonl>
//   highs-oracle basiclu-synth <dense70|rand120> <out.jsonl>
//   highs-oracle random  <out.jsonl>
//   highs-oracle ipx-model  <file.mps> <name> <cfg> <out.jsonl>
//   highs-oracle ipx-linalg <file.mps> <name> <cfg> <out.jsonl>
//   highs-oracle ipx-trace  <file.mps> <name> <cfg> <out.jsonl>
//   highs-oracle ipx        <file.mps> <name> <out.jsonl>
//   highs-oracle w1-qp|w1-ranging|w1-iis|w1-icrash <file.mps> <name> <out.jsonl>
//   highs-oracle w3-cupdlp|w3-hipdlp <file.mps> <name> <out.jsonl>
//   highs-oracle w4-orderings <file.mps> <name> <metis|amd|rcm> <out.jsonl>
//   highs-oracle w4-hipo <file.mps> <name> <out.jsonl>
//   highs-oracle w5-simplex <file.mps> <name> <2|3> <threads> <out.jsonl>
//   highs-oracle w5-hipo <file.mps> <name> <tree|node|both>
//                        <choose|augmented|normaleq> <threads> <out.jsonl>
//   highs-oracle g1-<x> <model> <name> <scenario> <out.jsonl>
//     (x: user-scale, inf-cost, semi, set-solution, multi-objective,
//      presolve-api, standard-form, rays, ill-conditioning, modify, mipfile)
#include <cstdio>
#include <cstdlib>
#include <string>

#include "oracle.h"

static int usage() {
  fprintf(stderr,
          "usage: highs-oracle unit|factor|basiclu <mps> <name> <out>\n"
          "       highs-oracle simplex <mps> <name> <strategy> <out>\n"
          "       highs-oracle basiclu-synth <dense70|rand120> <out>\n"
          "       highs-oracle random <out>\n"
          "       highs-oracle ipx-model|ipx-linalg|ipx-trace <mps> <name> <cfg> <out>\n"
          "       highs-oracle ipx <mps> <name> <out>\n"
          "       highs-oracle w1-qp|w1-ranging|w1-iis|w1-icrash <mps> <name> <out>\n"
          "       highs-oracle w3-cupdlp|w3-hipdlp <mps> <name> <out>\n"
          "       highs-oracle w4-orderings <mps> <name> <metis|amd|rcm> <out>\n"
          "       highs-oracle w4-hipo <mps> <name> <out>\n"
          "       highs-oracle w5-simplex <mps> <name> <2|3> <threads> <out>\n"
          "       highs-oracle w5-hipo <mps> <name> <tree|node|both> <choose|augmented|normaleq> <threads> <out>\n"
          "       highs-oracle g1-<x> <model> <name> <scenario> <out>\n");
  return 64;
}

int main(int argc, char** argv) {
  if (argc < 2) return usage();
  const std::string cmd = argv[1];
  if (cmd == "random" && argc == 3) return runRandom(argv[2]);
  if (cmd == "unit" && argc == 5) return runUnit(argv[2], argv[3], argv[4]);
  if (cmd == "factor" && argc == 5) return runFactor(argv[2], argv[3], argv[4]);
  if (cmd == "basiclu" && argc == 5) return runBasiclu(argv[2], argv[3], argv[4]);
  if (cmd == "basiclu-synth" && argc == 4) return runBasicluSynth(argv[2], argv[3]);
  if (cmd == "ipx-model" && argc == 6) return runIpxModel(argv[2], argv[3], argv[4], argv[5]);
  if (cmd == "ipx-linalg" && argc == 6) return runIpxLinalg(argv[2], argv[3], argv[4], argv[5]);
  if (cmd == "ipx-trace" && argc == 6) return runIpxTrace(argv[2], argv[3], argv[4], argv[5]);
  if (cmd == "ipx" && argc == 5) return runIpx(argv[2], argv[3], argv[4]);
  if (cmd == "w1-qp" && argc == 5) return runW1Qp(argv[2], argv[3], argv[4]);
  if (cmd == "w1-ranging" && argc == 5) return runW1Ranging(argv[2], argv[3], argv[4]);
  if (cmd == "w1-iis" && argc == 5) return runW1Iis(argv[2], argv[3], argv[4]);
  if (cmd == "w1-icrash" && argc == 5) return runW1Icrash(argv[2], argv[3], argv[4]);
  if (cmd == "w3-cupdlp" && argc == 5) return runW3Cupdlp(argv[2], argv[3], argv[4]);
  if (cmd == "w3-hipdlp" && argc == 5) return runW3Hipdlp(argv[2], argv[3], argv[4]);
  if (cmd == "w4-orderings" && argc == 6)
    return runW4Orderings(argv[2], argv[3], argv[4], argv[5]);
  if (cmd == "w4-hipo" && argc == 5) return runW4Hipo(argv[2], argv[3], argv[4]);
  if (cmd == "w5-simplex" && argc == 7)
    return runW5Simplex(argv[2], argv[3], argv[4], argv[5], argv[6]);
  if (cmd == "w5-hipo" && argc == 8)
    return runW5Hipo(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
  if (cmd.compare(0, 3, "g1-") == 0 && argc == 6) {
    const int r = runG1(cmd, argv[2], argv[3], argv[4], argv[5]);
    if (r >= 0) return r;
  }
  if (cmd == "simplex" && argc == 6)
    return runSimplex(argv[2], argv[3], atoi(argv[4]), argv[5]);
  return usage();
}
