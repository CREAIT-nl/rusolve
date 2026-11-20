// Generator of mt19937.out: libc++ std::mt19937 + uniform_real_distribution<double>(0,1),
// interleaved as feasibilityjump.hh draws them (seeds as the int `seed` it receives).
// clang++ -std=c++11 -O2 -ffp-contract=off mt19937.cpp -o mt19937 && ./mt19937
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
int main() {
  const int seeds[] = {0, 1, -1, 12345, 5489};
  for (int seed : seeds) {
    std::mt19937 rng(seed);
    uint64_t h = 0xcbf29ce484222325ull;
    for (int i = 0; i < 10000; ++i) {
      uint64_t v;
      if (i % 3 == 0) {
        double d = std::uniform_real_distribution<double>(0., 1.)(rng);
        std::memcpy(&v, &d, 8);
      } else {
        v = rng();
      }
      if (i < 4) printf("%d %d %016llx\n", seed, i, (unsigned long long)v);
      for (int b = 0; b < 8; ++b) { h ^= (v >> (8 * b)) & 0xff; h *= 0x100000001b3ull; }
    }
    printf("%d fnv %llu\n", seed, (unsigned long long)h);
  }
  return 0;
}
