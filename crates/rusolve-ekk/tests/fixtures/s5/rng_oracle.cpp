// RNG oracle for the mip-heuristics FPR / Local-MIP port.
//
// Prints the raw xoshiro256++ stream of mip-heuristics rng.h (commit 36f7d9d)
// and the outputs of every libstdc++ distribution / algorithm call shape the
// Rust port re-expresses, so the RNG test can compare bit for bit.
//
// Build (GCC 13.3.0, libstdc++ 13.3):
//   g++ -std=c++20 -O2 -ffp-contract=off -I <mip-heuristics>/src rng_oracle.cpp -o rng_oracle
//   ./rng_oracle > rng.txt
//
// Format: one line per shape, `<shape> seed=<s> <params> : <values>`.
// Every shape uses a fresh Rng(seed). Doubles are printed as the 16 hex digits
// of their IEEE-754 bits; integers in decimal.

#include "rng.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

namespace {

constexpr int kDraws = 64;

unsigned long long bits(double x) {
    std::uint64_t u;
    std::memcpy(&u, &x, sizeof u);
    return static_cast<unsigned long long>(u);
}

void raw(std::uint64_t seed) {
    Rng rng(seed);
    std::printf("raw seed=%llu :", static_cast<unsigned long long>(seed));
    for (int k = 0; k < kDraws; ++k) {
        std::printf(" %llu", static_cast<unsigned long long>(rng()));
    }
    std::printf("\n");
}

template <typename T>
void uniform_int(const char* tname, T a, T b) {
    Rng rng(42);
    std::uniform_int_distribution<T> dist(a, b);
    std::printf("uniform_int_%s seed=42 a=%lld b=%llu :", tname, static_cast<long long>(a),
                static_cast<unsigned long long>(b));
    for (int k = 0; k < kDraws; ++k) {
        const T v = dist(rng);
        if constexpr (std::is_signed_v<T>) {
            std::printf(" %lld", static_cast<long long>(v));
        } else {
            std::printf(" %llu", static_cast<unsigned long long>(v));
        }
    }
    std::printf("\n");
}

// Signed upper bounds are printed through %lld so they read naturally.
template <typename T>
void uniform_int_s(const char* tname, T a, T b) {
    Rng rng(42);
    std::uniform_int_distribution<T> dist(a, b);
    std::printf("uniform_int_%s seed=42 a=%lld b=%lld :", tname, static_cast<long long>(a),
                static_cast<long long>(b));
    for (int k = 0; k < kDraws; ++k) {
        std::printf(" %lld", static_cast<long long>(dist(rng)));
    }
    std::printf("\n");
}

void uniform_real(double a, double b) {
    Rng rng(42);
    std::uniform_real_distribution<double> dist(a, b);
    std::printf("uniform_real seed=42 a=%016llx b=%016llx :", bits(a), bits(b));
    for (int k = 0; k < kDraws; ++k) {
        std::printf(" %016llx", bits(dist(rng)));
    }
    std::printf("\n");
}

// The fresh-distribution-per-call shape upstream uses
// (`std::uniform_real_distribution<double>(lo, hi)(rng)`).
void uniform_real_temp(double a, double b) {
    Rng rng(42);
    std::printf("uniform_real_temp seed=42 a=%016llx b=%016llx :", bits(a), bits(b));
    for (int k = 0; k < kDraws; ++k) {
        std::printf(" %016llx", bits(std::uniform_real_distribution<double>(a, b)(rng)));
    }
    std::printf("\n");
}

void canonical() {
    Rng rng(42);
    std::printf("generate_canonical seed=42 :");
    for (int k = 0; k < kDraws; ++k) {
        const double x =
            std::generate_canonical<double, std::numeric_limits<double>::digits>(rng);
        std::printf(" %016llx", bits(x));
    }
    std::printf("\n");
}

void bernoulli(double p) {
    Rng rng(42);
    std::bernoulli_distribution dist(p);
    std::printf("bernoulli seed=42 p=%016llx :", bits(p));
    for (int k = 0; k < kDraws; ++k) {
        std::printf(" %d", dist(rng) ? 1 : 0);
    }
    std::printf("\n");
}

void shuffle(int n) {
    Rng rng(42);
    std::vector<int> v(static_cast<size_t>(n));
    std::iota(v.begin(), v.end(), 0);
    std::shuffle(v.begin(), v.end(), rng);
    std::printf("shuffle seed=42 n=%d :", n);
    for (int x : v) {
        std::printf(" %d", x);
    }
    // One more raw draw shows how many draws the shuffle consumed.
    std::printf(" | next=%llu\n", static_cast<unsigned long long>(rng()));
}

// std::ranges::sort (libstdc++ introsort) on (id, key) pairs with many tied
// keys, ordered by key descending as fpr_var_order.cpp:187 does.
void sort_ties(int n, int nkeys) {
    Rng rng(42);
    std::vector<std::pair<int, int>> v;
    for (int i = 0; i < n; ++i) {
        v.emplace_back(i, static_cast<int>(rng() % static_cast<std::uint64_t>(nkeys)));
    }
    std::printf("sort seed=42 n=%d nkeys=%d keys=", n, nkeys);
    for (int i = 0; i < n; ++i) {
        std::printf("%s%d", i == 0 ? "" : ",", v[i].second);
    }
    std::ranges::sort(v, [](const auto& a, const auto& c) { return a.second > c.second; });
    std::printf(" :");
    for (const auto& [id, key] : v) {
        std::printf(" %d", id);
    }
    std::printf("\n");
}

// std::partial_sort with tied keys, as local_mip_search.cpp:475 sorts sampled
// constraints by weight descending (recorded for the Local-MIP search port).
void partial_sort_ties(int n, int keep, int nkeys) {
    Rng rng(42);
    std::vector<std::pair<int, std::uint64_t>> v;
    for (int i = 0; i < n; ++i) {
        v.emplace_back(i, 1 + rng() % static_cast<std::uint64_t>(nkeys));
    }
    std::printf("partial_sort seed=42 n=%d keep=%d nkeys=%d keys=", n, keep, nkeys);
    for (int i = 0; i < n; ++i) {
        std::printf("%s%llu", i == 0 ? "" : ",", static_cast<unsigned long long>(v[i].second));
    }
    std::partial_sort(v.begin(), v.begin() + keep, v.end(),
                      [](const auto& a, const auto& b) { return a.second > b.second; });
    std::printf(" :");
    for (const auto& [id, w] : v) {
        std::printf(" %d", id);
    }
    std::printf("\n");
}

}  // namespace

int main() {
    raw(0);
    raw(42);
    raw(1039);

    // Ranges 1, 2, 3, 2^31, 2^63 and full width, per type.
    uniform_int_s<std::int32_t>("i32", 0, 0);
    uniform_int_s<std::int32_t>("i32", 0, 1);
    uniform_int_s<std::int32_t>("i32", 0, 2);
    uniform_int_s<std::int32_t>("i32", -1, 1);
    uniform_int_s<std::int32_t>("i32", 0, std::numeric_limits<std::int32_t>::max());
    uniform_int_s<std::int32_t>("i32", std::numeric_limits<std::int32_t>::min(),
                                std::numeric_limits<std::int32_t>::max());
    uniform_int_s<std::int64_t>("i64", 0, 0);
    uniform_int_s<std::int64_t>("i64", 0, 1);
    uniform_int_s<std::int64_t>("i64", 1, 3);
    uniform_int_s<std::int64_t>("i64", 0, (std::int64_t{1} << 31) - 1);
    uniform_int_s<std::int64_t>("i64", 0, std::numeric_limits<std::int64_t>::max());
    uniform_int_s<std::int64_t>("i64", std::numeric_limits<std::int64_t>::min(),
                                std::numeric_limits<std::int64_t>::max());
    uniform_int<std::uint64_t>("u64", 0, 0);
    uniform_int<std::uint64_t>("u64", 1, 2);
    uniform_int<std::uint64_t>("u64", 0, 2);
    uniform_int<std::uint64_t>("u64", 0, (std::uint64_t{1} << 31) - 1);
    uniform_int<std::uint64_t>("u64", 0, (std::uint64_t{1} << 63) - 1);
    uniform_int<std::uint64_t>("u64", 0, std::numeric_limits<std::uint64_t>::max());

    uniform_real(0.0, 1.0);
    uniform_real(1e-15, 1.0);
    uniform_real(-3.5, 12.25);
    {
        const double range = 7.3;
        uniform_real(-0.1 * range, 0.1 * range);
    }
    uniform_real_temp(-3.5, 12.25);
    canonical();
    bernoulli(0.5);

    for (int n : {0, 1, 2, 3, 30, 31, 1000}) {
        shuffle(n);
    }

    sort_ties(20, 3);
    sort_ties(100, 4);
    partial_sort_ties(20, 12, 3);
    partial_sort_ties(100, 12, 4);
    return 0;
}
