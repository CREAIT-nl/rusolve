// Shared declarations for the highs-oracle subcommands.
#ifndef HIGHS_ORACLE_ORACLE_H_
#define HIGHS_ORACLE_ORACLE_H_

#include <string>
#include <vector>

#include "Highs.h"
#include "emit.h"

// Read an MPS file into a HighsLp exactly as `Highs::readModel` stores it
// (column-wise, unscaled). Integrality is dropped (LP relaxation).
HighsLp oracleLoadLp(const std::string& path);

// Common option set for every oracle Highs instance: silent, serial,
// deterministic. Callers add the per-fixture options on top.
void oracleBaseOptions(Highs& highs, int random_seed);

// The scaled LP and basic index HEkk holds when HEkk::solve returns, for a
// dual simplex run with the given iteration limit (-1 = none).
struct EkkSnapshot {
  HighsLp scaled_lp;  // HEkk::lp_ at exit (scaled, column-wise)
  std::vector<HighsInt> basic_index;
  std::vector<int8_t> nonbasic_flag;
  HighsInt iteration_count = 0;
  bool valid = false;
};
EkkSnapshot oracleSnapshot(const std::string& path, HighsInt iteration_limit);

int runUnit(const std::string& mps, const std::string& name, const std::string& out);
int runFactor(const std::string& mps, const std::string& name, const std::string& out);
int runBasiclu(const std::string& mps, const std::string& name, const std::string& out);
int runBasicluSynth(const std::string& name, const std::string& out);
int runSimplex(const std::string& mps, const std::string& name, int strategy,
               const std::string& out);
// runSimplex at threads=N with parallel=on (N > 0; 0 = runSimplex itself).
int runSimplexThreads(const std::string& mps, const std::string& name, int strategy,
                      int threads, const std::string& out);
int runRandom(const std::string& out);
// IPX (ipx.cpp): model after Load, the linear-algebra components, and the
// end-to-end IPM trace of Highs with solver=ipx (needs the IPX hooks patch).
int runIpxModel(const std::string& mps, const std::string& name, const std::string& cfg,
                const std::string& out);
int runIpxLinalg(const std::string& mps, const std::string& name, const std::string& cfg,
                 const std::string& out);
int runIpxTrace(const std::string& mps, const std::string& name, const std::string& cfg,
                const std::string& out);
int runIpx(const std::string& mps, const std::string& name, const std::string& out);
// W1 (w1.cpp): QUASS (needs the W1 hooks patch), getRanging, getIis, iCrash
// through the public Highs API (FORMAT.md "W1 records").
int runW1Qp(const std::string& mps, const std::string& name, const std::string& out);
int runW1Ranging(const std::string& mps, const std::string& name, const std::string& out);
int runW1Iis(const std::string& mps, const std::string& name, const std::string& out);
int runW1Icrash(const std::string& mps, const std::string& name, const std::string& out);
// W3 (w3.cpp): cuPDLP-C (solver=pdlp) and HiPDLP (solver=hipdlp) through
// Highs::run, traced by the W3 hooks patch (FORMAT.md "W3 records").
int runW3Cupdlp(const std::string& mps, const std::string& name, const std::string& out);
int runW3Hipdlp(const std::string& mps, const std::string& name, const std::string& out);
// W4 (w4.cpp): HiPO (solver=hipo) through Highs::run in the HIPO=ON build
// (refblas), traced by the W4 hooks patch (FORMAT.md "W4 records").
int runW4Orderings(const std::string& mps, const std::string& name, const std::string& ordering,
                   const std::string& out);
int runW4Hipo(const std::string& mps, const std::string& name, const std::string& out);
// The W4 run with the given overrides (k=v,k=v); threads > 0 adds the W5
// "threads" / "parallel" header keys (w5-hipo), 0 writes the W4 header.
int runW4WithOverrides(const std::string& mps, const std::string& name,
                       const std::string& path, const std::string& overrides,
                       int threads);
// W5 (w5.cpp): the simplex trace with simplex_strategy 2 (SIP) / 3 (PAMI)
// and the HiPO trace, both at threads=N with parallel=on (FORMAT.md "W5
// records"); w5-hipo needs the HIPO=ON build and the W5 hooks patch.
int runW5Simplex(const std::string& mps, const std::string& name, const std::string& strategy,
                 const std::string& threads, const std::string& out);
int runW5Hipo(const std::string& mps, const std::string& name, const std::string& type,
              const std::string& system, const std::string& threads, const std::string& out);
// G1 (g1.cpp): the Highs API around run through the public Highs API, no
// hooks (FORMAT.md "G1 records"). cmd is the g1-* command name; returns -1
// for an unknown one.
int runG1(const std::string& cmd, const std::string& model, const std::string& name,
          const std::string& scenario, const std::string& out);
// The G1 model reader and "lp" + "integrality" writer, for G1 commands
// outside g1.cpp.
HighsLp oracleG1Read(const std::string& path);
void oracleG1AddLp(Obj& o, const HighsLp& lp);
// g1-clique-normalise (clqn.cpp): HPresolve::normaliseCliqueRows on the
// model, the presolve state after it (HiGHS 1.16 only; 64 on 1.15.1).
int runClqn(const std::string& model, const std::string& name, const std::string& scenario,
            const std::string& out);

#endif
