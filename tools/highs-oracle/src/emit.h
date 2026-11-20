// JSONL emitter for the rusolve-ekk golden fixtures. See ../FORMAT.md.
//
// Every f64 is written as the 16 lowercase hex digits of its IEEE-754 bit
// pattern (so -0.0, NaN payloads, inf and subnormals survive exactly). Every
// integer (HighsInt, int8_t, bool, enum) is a plain JSON number.
#ifndef HIGHS_ORACLE_EMIT_H_
#define HIGHS_ORACLE_EMIT_H_

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "util/HVector.h"
#include "util/HighsSparseMatrix.h"

std::string hexf(double v);
std::string hex64(uint64_t v);

// FNV-1a 64 over little-endian bytes: f64 as its 8 bit-pattern bytes,
// HighsInt as 4 bytes (HighsInt is int32 in this build), int8 as 1 byte.
uint64_t fnvF(const double* p, size_t n);
uint64_t fnvI(const HighsInt* p, size_t n);
uint64_t fnvI8(const int8_t* p, size_t n);
inline uint64_t fnvF(const std::vector<double>& v) { return fnvF(v.data(), v.size()); }
inline uint64_t fnvI(const std::vector<HighsInt>& v) { return fnvI(v.data(), v.size()); }
inline uint64_t fnvI8(const std::vector<int8_t>& v) { return fnvI8(v.data(), v.size()); }

class Obj {
 public:
  Obj() {}
  explicit Obj(const char* type) { s("t", type); }
  Obj& i(const char* k, long long v);
  Obj& b(const char* k, bool v) { return i(k, v ? 1 : 0); }
  Obj& f(const char* k, double v);
  Obj& h(const char* k, uint64_t v);  // hash, as 16 hex digits
  Obj& s(const char* k, const std::string& v);
  Obj& vi(const char* k, const HighsInt* p, size_t n);
  Obj& vi(const char* k, const std::vector<HighsInt>& v) { return vi(k, v.data(), v.size()); }
  Obj& vi8(const char* k, const int8_t* p, size_t n);
  Obj& vi8(const char* k, const std::vector<int8_t>& v) { return vi8(k, v.data(), v.size()); }
  Obj& vf(const char* k, const double* p, size_t n);
  Obj& vf(const char* k, const std::vector<double>& v) { return vf(k, v.data(), v.size()); }
  Obj& raw(const char* k, const std::string& json);
  Obj& o(const char* k, const Obj& sub) { return raw(k, sub.str()); }
  // HVector: size, count, index[0..count) (exact order), synthetic_tick, and
  // every array position whose bit pattern is non-zero (nz_i ascending,
  // nz_v), so -0.0 and kHighsZero (1e-50) entries are preserved. With
  // with_pack also packFlag/packCount/packIndex/packValue.
  Obj& hv(const char* k, const HVector& v, bool with_pack = false);
  // HighsSparseMatrix: format, num_col, num_row, start, p_end (if partitioned),
  // index, value.
  Obj& mat(const char* k, const HighsSparseMatrix& m);
  std::string str() const { return "{" + body_ + "}"; }

 private:
  void key(const char* k);
  std::string body_;
  bool first_ = true;
};

class Out {
 public:
  explicit Out(const std::string& path);
  ~Out();
  void put(const Obj& o);
  // Lines already in JSONL form (e.g. the IPX hook records), copied verbatim.
  void raw_lines(const std::string& lines);
  void flush();
  bool ok() const { return f_ != nullptr; }

 private:
  FILE* f_;
};

#endif
