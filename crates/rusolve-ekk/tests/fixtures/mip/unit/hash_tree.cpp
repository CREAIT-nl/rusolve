// Generator of hash_tree.out: the real HighsHashTree (HiGHS 1.15.1, util/HighsHashTree.h).
// R=target/highs-oracle; clang++ -std=c++11 -O1 -I $R/highs-src/highs -I $R/highs-src/extern \
//   -I $R/build/highs-build hash_tree.cpp $R/build/lib/libhighs.a -o hash_tree && ./hash_tree
// The Rust side is util/highs_hash_tree.rs::tests::differential_vs_cpp_highs_hash_tree.
#include <cstdio>
#include <cstdint>
#include <vector>
#include "util/HighsHashTree.h"

static uint64_t st;
static uint64_t nxt() { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }

struct Fnv { uint64_t h = 0xcbf29ce484222325ull; void add(uint64_t v) { for (int i = 0; i < 8; ++i) { h ^= (v >> (8*i)) & 0xff; h *= 0x100000001b3ull; } } };

int main() {
  // Scenario A: random insert/erase with grow and shrink phases.
  st = 0x9e3779b97f4a7c15ull;
  HighsHashTree<HighsInt, HighsInt> t;
  for (int step = 0; step < 60000; ++step) {
    uint64_t r = nxt();
    HighsInt key = (HighsInt)(nxt() % 3000);
    int phase = (step / 10000) % 2;  // 0: grow-biased, 1: shrink-biased
    bool ins = phase == 0 ? (r % 4 != 0) : (r % 4 == 0);
    if (ins) t.insert(key, step); else t.erase(key);
    if (step % 2000 == 1999) {
      Fnv f; int n = 0;
      t.for_each([&](HighsInt k, HighsInt v) { f.add((uint64_t)(uint32_t)k); f.add((uint64_t)(uint32_t)v); ++n; });
      printf("A %d %d %llu\n", step, n, (unsigned long long)f.h);
    }
  }
  // Scenario B: small void tree full order.
  st = 12345;
  HighsHashTree<HighsInt> s;
  for (int i = 0; i < 120; ++i) s.insert((HighsInt)(nxt() % 100000));
  for (int i = 0; i < 40; ++i) s.erase((HighsInt)(nxt() % 100000));
  printf("B");
  s.for_each([&](HighsInt k) { printf(" %d", k); });
  printf("\n");
  // Scenario C: find_common between random trees of varying sizes.
  st = 777;
  Fnv fc;
  for (int trial = 0; trial < 400; ++trial) {
    HighsHashTree<HighsInt, HighsInt> t1, t2;
    int n1 = 1 + nxt() % 200, n2 = 1 + nxt() % 200;
    for (int i = 0; i < n1; ++i) t1.insert((HighsInt)(nxt() % 5000), i);
    for (int i = 0; i < n2; ++i) t2.insert((HighsInt)(nxt() % 5000), -i);
    const HighsHashTableEntry<HighsInt, HighsInt>* e = t1.find_common(t2);
    long long out = e ? ((long long)e->key() * 1000 + e->value()) : -1;
    fc.add((uint64_t)out);
    if (trial < 10) printf("C %d %lld\n", trial, out);
  }
  printf("Cfnv %llu\n", (unsigned long long)fc.h);
  // Scenario D: fill then drain in pseudo-random order (exercises branch merges).
  st = 4242;
  HighsHashTree<HighsInt, HighsInt> d;
  std::vector<HighsInt> keys;
  for (int i = 0; i < 1500; ++i) { HighsInt k = (HighsInt)(nxt() % 1000000); if (d.insert(k, i)) keys.push_back(k); }
  for (size_t i = keys.size(); i > 1; --i) { size_t j = nxt() % i; std::swap(keys[i-1], keys[j]); }
  for (size_t i = 0; i < keys.size(); ++i) {
    d.erase(keys[i]);
    size_t left = keys.size() - i - 1;
    if (left % 50 == 0 || left < 60) {
      Fnv f; int n = 0;
      d.for_each([&](HighsInt k, HighsInt v) { f.add((uint64_t)(uint32_t)k); f.add((uint64_t)(uint32_t)v); ++n; });
      printf("D %zu %d %llu\n", left, n, (unsigned long long)f.h);
    }
  }
  printf("Dempty %d\n", (int)d.empty());
  return 0;
}
