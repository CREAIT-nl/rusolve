// W3 oracle commands: cuPDLP-C (w3-cupdlp, solver=pdlp) and HiPDLP (w3-hipdlp,
// solver=hipdlp) through the public Highs API on a fresh Highs with the W3
// base options. Records: FORMAT.md "W3 records".
//
//   highs-oracle w3-cupdlp <file.mps> <name> <out.jsonl>
//   highs-oracle w3-hipdlp <file.mps> <name> <out.jsonl>
//
// Both need the fourth patch, highs-oracle-w3-hooks.patch
// (pdlp/HighsOracleW3Dump.h): the command writes the w3_header line, points
// HIGHS_ORACLE_PDLP_DUMP at the fixture itself (the hooks append), runs, and
// appends w3_result, so the pdlp_* / hipdlp_* hook lines land in place.
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "Highs.h"
#include "oracle.h"
#include "pdlp/HighsOracleW3Dump.h"

namespace {

// The W3 base options (tests/port/w3.rs w3_options): output off, serial,
// presolve off, random_seed=0; the solver comes from the overrides.
void w3BaseOptions(Highs& highs) {
  highs.setOptionValue("output_flag", false);
  highs.setOptionValue("threads", 1);
  highs.setOptionValue("parallel", "off");
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("random_seed", 0);
}

bool w3Overrides(Highs& highs, const std::string& cfg) {
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
void w3Read(const std::string& mps, HighsLp& lp) {
  Highs reader;
  w3BaseOptions(reader);
  if (reader.readModel(mps) == HighsStatus::kError) {
    fprintf(stderr, "highs-oracle: cannot read %s\n", mps.c_str());
    exit(2);
  }
  lp = reader.getLp();
  lp.integrality_.clear();
}

// HighsLp, encoded as the "lp" object of FORMAT.md (unscaled).
Obj w3LpObj(const HighsLp& lp) {
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

Obj w3Header(const std::string& name, const std::string& overrides, const HighsLp& lp) {
  Obj h("w3_header");
  h.s("name", name).s("highs_version", highsVersion()).s("highs_githash", highsGithash())
      .s("overrides", overrides).o("lp", w3LpObj(lp));
  return h;
}

Obj w3Result(HighsStatus run_status, const Highs& highs) {
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& sol = highs.getSolution();
  Obj r("w3_result");
  r.i("run_status", (int)run_status)
      .i("model_status", (int)highs.getModelStatus())
      .f("objective_function_value", info.objective_function_value)
      .i("pdlp_iteration_count", info.pdlp_iteration_count)
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

bool w3Append(const std::string& path, const Obj& o) {
  std::FILE* f = std::fopen(path.c_str(), "ab");
  if (f == nullptr) return false;
  const std::string line = o.str() + "\n";
  const bool ok = std::fwrite(line.data(), 1, line.size(), f) == line.size();
  return std::fclose(f) == 0 && ok;
}

// w3-cupdlp / w3-hipdlp: header, hooked run, result.
int runW3(const std::string& mps, const std::string& name, const std::string& path,
          const std::string& overrides) {
  HighsLp lp;
  w3Read(mps, lp);
  {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
      fprintf(stderr, "highs-oracle: cannot open %s\n", path.c_str());
      return 2;
    }
    const std::string h = w3Header(name, overrides, lp).str() + "\n";
    std::fwrite(h.data(), 1, h.size(), f);
    if (std::fclose(f) != 0) return 2;
  }
  Highs highs;
  w3BaseOptions(highs);
  if (!w3Overrides(highs, overrides)) return 64;
  if (highs.passModel(lp) != HighsStatus::kOk) {
    fprintf(stderr, "highs-oracle: passModel failed for %s\n", name.c_str());
    return 2;
  }
  highsOracleW3Close();
  *highsOracleW3Call() = 0;
  setenv("HIGHS_ORACLE_PDLP_DUMP", path.c_str(), 1);
  const HighsStatus run_status = highs.run();
  const long long hook_lines = *highsOracleW3Call();
  highsOracleW3Close();  // flushes and closes the hook file
  unsetenv("HIGHS_ORACLE_PDLP_DUMP");
  // Every PDLP run reaches its end hook, so an empty trace means the hooks
  // are not compiled in (or the solver was not reached).
  if (hook_lines == 0) {
    fprintf(stderr, "highs-oracle: no PDLP hook line written for %s\n", name.c_str());
    return 3;
  }
  if (!w3Append(path, w3Result(run_status, highs))) return 2;
  return 0;
}

}  // namespace

int runW3Cupdlp(const std::string& mps, const std::string& name, const std::string& path) {
  return runW3(mps, name, path, "solver=pdlp,pdlp_iteration_limit=2000");
}

int runW3Hipdlp(const std::string& mps, const std::string& name, const std::string& path) {
  return runW3(mps, name, path, "solver=hipdlp,pdlp_iteration_limit=2000");
}
