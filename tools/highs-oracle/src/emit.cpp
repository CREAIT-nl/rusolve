#include "emit.h"

#include <cinttypes>
#include <cstdlib>
#include <cstring>

std::string hex64(uint64_t v) {
  char buf[17];
  snprintf(buf, sizeof buf, "%016" PRIx64, v);
  return std::string(buf, 16);
}

std::string hexf(double v) {
  uint64_t u;
  memcpy(&u, &v, 8);
  return hex64(u);
}

static const uint64_t kFnvOffset = 0xcbf29ce484222325ULL;
static const uint64_t kFnvPrime = 0x100000001b3ULL;

static uint64_t fnvBytes(uint64_t h, const unsigned char* p, size_t n) {
  for (size_t i = 0; i < n; i++) {
    h ^= p[i];
    h *= kFnvPrime;
  }
  return h;
}

uint64_t fnvF(const double* p, size_t n) {
  uint64_t h = kFnvOffset;
  for (size_t i = 0; i < n; i++) {
    uint64_t u;
    memcpy(&u, &p[i], 8);
    unsigned char b[8];
    for (int k = 0; k < 8; k++) b[k] = (unsigned char)(u >> (8 * k));
    h = fnvBytes(h, b, 8);
  }
  return h;
}

uint64_t fnvI(const HighsInt* p, size_t n) {
  static_assert(sizeof(HighsInt) == 4, "oracle assumes 32-bit HighsInt");
  uint64_t h = kFnvOffset;
  for (size_t i = 0; i < n; i++) {
    uint32_t u = (uint32_t)p[i];
    unsigned char b[4];
    for (int k = 0; k < 4; k++) b[k] = (unsigned char)(u >> (8 * k));
    h = fnvBytes(h, b, 4);
  }
  return h;
}

uint64_t fnvI8(const int8_t* p, size_t n) {
  return fnvBytes(kFnvOffset, reinterpret_cast<const unsigned char*>(p), n);
}

void Obj::key(const char* k) {
  if (!first_) body_ += ',';
  first_ = false;
  body_ += '"';
  body_ += k;
  body_ += "\":";
}

Obj& Obj::i(const char* k, long long v) {
  key(k);
  body_ += std::to_string(v);
  return *this;
}

Obj& Obj::f(const char* k, double v) {
  key(k);
  body_ += '"';
  body_ += hexf(v);
  body_ += '"';
  return *this;
}

Obj& Obj::h(const char* k, uint64_t v) {
  key(k);
  body_ += '"';
  body_ += hex64(v);
  body_ += '"';
  return *this;
}

Obj& Obj::s(const char* k, const std::string& v) {
  key(k);
  body_ += '"';
  for (char c : v) {
    if (c == '"' || c == '\\') {
      body_ += '\\';
      body_ += c;
    } else if ((unsigned char)c < 0x20) {
      char buf[8];
      snprintf(buf, sizeof buf, "\\u%04x", (unsigned)c);
      body_ += buf;
    } else {
      body_ += c;
    }
  }
  body_ += '"';
  return *this;
}

Obj& Obj::vi(const char* k, const HighsInt* p, size_t n) {
  key(k);
  body_ += '[';
  for (size_t j = 0; j < n; j++) {
    if (j) body_ += ',';
    body_ += std::to_string((long long)p[j]);
  }
  body_ += ']';
  return *this;
}

Obj& Obj::vi8(const char* k, const int8_t* p, size_t n) {
  key(k);
  body_ += '[';
  for (size_t j = 0; j < n; j++) {
    if (j) body_ += ',';
    body_ += std::to_string((int)p[j]);
  }
  body_ += ']';
  return *this;
}

Obj& Obj::vf(const char* k, const double* p, size_t n) {
  key(k);
  body_ += '[';
  for (size_t j = 0; j < n; j++) {
    if (j) body_ += ',';
    body_ += '"';
    body_ += hexf(p[j]);
    body_ += '"';
  }
  body_ += ']';
  return *this;
}

Obj& Obj::raw(const char* k, const std::string& json) {
  key(k);
  body_ += json;
  return *this;
}

Obj& Obj::hv(const char* k, const HVector& v, bool with_pack) {
  Obj o;
  o.i("size", v.size).i("count", v.count);
  if (v.count >= 0 && (size_t)v.count <= v.index.size())
    o.vi("index", v.index.data(), v.count);
  o.f("synthetic_tick", v.synthetic_tick);
  std::vector<HighsInt> nz_i;
  std::vector<double> nz_v;
  for (size_t j = 0; j < v.array.size(); j++) {
    uint64_t u;
    memcpy(&u, &v.array[j], 8);
    if (u != 0) {
      nz_i.push_back((HighsInt)j);
      nz_v.push_back(v.array[j]);
    }
  }
  o.vi("nz_i", nz_i).vf("nz_v", nz_v);
  if (with_pack) {
    o.b("packFlag", v.packFlag).i("packCount", v.packCount);
    o.vi("packIndex", v.packIndex.data(), v.packCount);
    o.vf("packValue", v.packValue.data(), v.packCount);
  }
  return raw(k, o.str());
}

Obj& Obj::mat(const char* k, const HighsSparseMatrix& m) {
  Obj o;
  o.i("format", (int)m.format_).i("num_col", m.num_col_).i("num_row", m.num_row_);
  o.vi("start", m.start_);
  if (m.format_ == MatrixFormat::kRowwisePartitioned) o.vi("p_end", m.p_end_);
  o.vi("index", m.index_).vf("value", m.value_);
  return raw(k, o.str());
}

Out::Out(const std::string& path) : f_(fopen(path.c_str(), "wb")) {
  if (!f_) {
    fprintf(stderr, "highs-oracle: cannot open %s\n", path.c_str());
    exit(2);
  }
}

Out::~Out() {
  if (f_) fclose(f_);
}

void Out::put(const Obj& o) {
  std::string s = o.str();
  s += '\n';
  fwrite(s.data(), 1, s.size(), f_);
}

void Out::raw_lines(const std::string& lines) {
  fwrite(lines.data(), 1, lines.size(), f_);
}

void Out::flush() { fflush(f_); }
