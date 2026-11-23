// W4 oracle commands: HiPO through the public Highs API on a fresh Highs with
// the W4 base options, one forced fill-reducing ordering (w4-orderings) or
// HiPO's own choices (w4-hipo). Records: FORMAT.md "W4 records".
//
//   highs-oracle w4-orderings <file.mps> <name> <metis|amd|rcm> <out.jsonl>
//   highs-oracle w4-hipo      <file.mps> <name> <out.jsonl>
//
// Both need the HIPO=ON oracle configuration (CMakeLists.txt ORACLE_HIPO,
// gen_fixtures.sh build-hipo: HiPO linked against refblas) and the fifth
// patch, highs-oracle-w4-hooks.patch (ipm/hipo/HighsOracleW4Dump.h): the
// command checks that HiPO is available and runs on refblas, writes the
// w4_header line, points HIGHS_ORACLE_HIPO_DUMP at the fixture itself (the
// hooks append), runs, and appends w4_result, so the hipo_* hook lines land
// in place.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "Highs.h"
#include "HighsExternalApi.h"
#include "ipm/hipo/HighsOracleW4Dump.h"
#include "oracle.h"

namespace {

// HighsExtras::blas::getInfo()->provider of the pinned BLAS
// (CMakeLists.txt: HIGHS_BLAS_VENDOR).
const char* const kW4BlasVendor = "rusolve refblas";

// The W4 base options (tests/port/w4.rs w4_options): output off, serial,
// presolve off, random_seed=0; the solver comes from the overrides.
void w4BaseOptions(Highs& highs) {
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("threads", 1);
  highs.setOptionValue("parallel", "off");
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("random_seed", 0);
}

bool w4Overrides(Highs& highs, const std::string& cfg) {
  if (cfg == "default") return true;
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

// The LP as Highs::readModel stores it; integrality dropped.
void w4Read(const std::string& mps, HighsLp& lp) {
  Highs reader;
  w4BaseOptions(reader);
  if (reader.readModel(mps) == HighsStatus::kError) {
    fprintf(stderr, "highs-oracle: cannot read %s\n", mps.c_str());
    exit(2);
  }
  lp = reader.getLp();
  lp.integrality_.clear();
}

// HighsLp, encoded as the "lp" object of FORMAT.md (unscaled).
Obj w4LpObj(const HighsLp& lp) {
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

// threads > 0 (w5-hipo): also "threads" and "parallel" (FORMAT.md "W5
// records"); the W4 header has neither.
Obj w4Header(const std::string& name, const std::string& overrides, const HighsLp& lp,
             int threads) {
  Obj h("w4_header");
  h.s("name", name).s("highs_version", highsVersion()).s("highs_githash", highsGithash())
      .s("overrides", overrides);
  if (threads > 0) h.i("threads", threads).s("parallel", "on");
  h.o("lp", w4LpObj(lp));
  return h;
}

Obj w4Result(HighsStatus run_status, const Highs& highs) {
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& sol = highs.getSolution();
  Obj r("w4_result");
  r.i("run_status", (int)run_status)
      .i("model_status", (int)highs.getModelStatus())
      .f("objective_function_value", info.objective_function_value)
      .i("ipm_iteration_count", info.ipm_iteration_count)
      .i("crossover_iteration_count", info.crossover_iteration_count)
      .i("primal_solution_status", info.primal_solution_status)
      .i("dual_solution_status", info.dual_solution_status)
      .i("num_primal_infeasibilities", info.num_primal_infeasibilities)
      .f("max_primal_infeasibility", info.max_primal_infeasibility)
      .f("sum_primal_infeasibilities", info.sum_primal_infeasibilities)
      .i("num_dual_infeasibilities", info.num_dual_infeasibilities)
      .f("max_dual_infeasibility", info.max_dual_infeasibility)
      .f("sum_dual_infeasibilities", info.sum_dual_infeasibilities)
      .b("basis_valid", highs.getBasis().valid)
      .vf("col_value", sol.col_value).vf("col_dual", sol.col_dual)
      .vf("row_value", sol.row_value).vf("row_dual", sol.row_dual);
  return r;
}

bool w4Append(const std::string& path, const Obj& o) {
  std::FILE* f = std::fopen(path.c_str(), "ab");
  if (f == nullptr) return false;
  const std::string line = o.str() + "\n";
  const bool ok = std::fwrite(line.data(), 1, line.size(), f) == line.size();
  return std::fclose(f) == 0 && ok;
}

// A W4 fixture is only meaningful from the HIPO=ON build on refblas: in the
// default build HiGHS would silently fall back to IPX, and another BLAS would
// give different (unportable) bits.
bool w4CheckHipo() {
  if (!HighsExternalApi::isAvailable<HighsExtras::hipo>()) {
    fprintf(stderr,
            "highs-oracle: HiPO is not available in this build (use the "
            "ORACLE_HIPO=ON build, gen_fixtures.sh build-hipo)\n");
    return false;
  }
  const char* provider = HighsExtras::blas::getInfo()->provider;
  if (provider == nullptr || std::strcmp(provider, kW4BlasVendor) != 0) {
    fprintf(stderr, "highs-oracle: HiPO BLAS is '%s', expected '%s'\n",
            provider == nullptr ? "(null)" : provider, kW4BlasVendor);
    return false;
  }
  return true;
}

}  // namespace

// w4-orderings / w4-hipo / w5-hipo: header, hooked run, result.
int runW4WithOverrides(const std::string& mps, const std::string& name,
                       const std::string& path, const std::string& overrides,
                       int threads) {
  if (!w4CheckHipo()) return 3;
  HighsLp lp;
  w4Read(mps, lp);
  {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
      fprintf(stderr, "highs-oracle: cannot open %s\n", path.c_str());
      return 2;
    }
    const std::string h = w4Header(name, overrides, lp, threads).str() + "\n";
    std::fwrite(h.data(), 1, h.size(), f);
    if (std::fclose(f) != 0) return 2;
  }
  Highs highs;
  w4BaseOptions(highs);
  if (!w4Overrides(highs, overrides)) return 64;
  if (highs.passModel(lp) != HighsStatus::kOk) {
    fprintf(stderr, "highs-oracle: passModel failed for %s\n", name.c_str());
    return 2;
  }
  highsOracleW4Close();
  *highsOracleW4Call() = 0;
  setenv("HIGHS_ORACLE_HIPO_DUMP", path.c_str(), 1);
  const HighsStatus run_status = highs.run();
  const long long hook_lines = *highsOracleW4Call();
  highsOracleW4Close();  // flushes and closes the hook file
  unsetenv("HIGHS_ORACLE_HIPO_DUMP");
  // Every HiPO run reaches its hipo_end hook, so an empty trace means the
  // hooks are not compiled in or HiGHS did not take the HiPO path.
  if (hook_lines == 0) {
    fprintf(stderr, "highs-oracle: no HiPO hook line written for %s\n", name.c_str());
    return 3;
  }
  if (!w4Append(path, w4Result(run_status, highs))) return 2;
  return 0;
}

int runW4Orderings(const std::string& mps, const std::string& name,
                   const std::string& ordering, const std::string& path) {
  if (ordering != "metis" && ordering != "amd" && ordering != "rcm") {
    fprintf(stderr, "highs-oracle: w4-orderings: ordering must be metis, amd or rcm\n");
    return 64;
  }
  return runW4WithOverrides(mps, name, path, "solver=hipo,hipo_ordering=" + ordering, 0);
}

int runW4Hipo(const std::string& mps, const std::string& name, const std::string& path) {
  return runW4WithOverrides(mps, name, path, "solver=hipo", 0);
}
