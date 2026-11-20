// Generator of cpp_std.out: libc++ heap / partial_sort / partition /
// priority_queue order on short inputs with ties, for util/cpp_std.rs.
// Inputs: every key array in {0,1,2}^n for n = 0..5; element = (key, id = position);
// the comparator orders on key only, so the id sequence exposes the tie order.
// clang++ -std=c++17 -O2 cpp_std.cpp -o cpp_std && ./cpp_std > cpp_std.out
#include <algorithm>
#include <cstdio>
#include <queue>
#include <utility>
#include <vector>
typedef std::pair<int, int> E;  // (key, id)
static bool lt(const E& a, const E& b) { return a.first < b.first; }
static void out(const char* op, int n, int code, int extra, const std::vector<E>& v) {
  std::printf("%s %d %d %d", op, n, code, extra);
  for (const E& e : v) std::printf(" %d", e.second);
  std::printf("\n");
}
int main() {
  for (int n = 0; n <= 5; ++n) {
    int total = 1;
    for (int i = 0; i < n; ++i) total *= 3;
    for (int code = 0; code < total; ++code) {
      std::vector<E> base(n);
      int c = code;
      for (int i = 0; i < n; ++i) { base[i] = E(c % 3, i); c /= 3; }

      std::vector<E> v = base;
      std::make_heap(v.begin(), v.end(), lt);
      out("make", n, code, 0, v);

      v = base;
      for (int i = 1; i <= n; ++i) std::push_heap(v.begin(), v.begin() + i, lt);
      out("push", n, code, 0, v);

      v = base;
      std::make_heap(v.begin(), v.end(), lt);
      for (int k = n; k >= 1; --k) std::pop_heap(v.begin(), v.begin() + k, lt);
      out("pop", n, code, 0, v);

      v = base;
      std::make_heap(v.begin(), v.end(), lt);
      std::sort_heap(v.begin(), v.end(), lt);
      out("sortheap", n, code, 0, v);

      for (int m = 0; m <= n; ++m) {
        v = base;
        std::partial_sort(v.begin(), v.begin() + m, v.end(), lt);
        out("psort", n, code, m, v);
      }

      for (int p = 0; p < 3; ++p) {
        v = base;
        auto it = std::partition(v.begin(), v.end(), [&](const E& e) { return e.first == p; });
        out("part", n, code, p * 100 + int(it - v.begin()), v);
      }

      // priority_queue ordered on .first only (as in HighsGFkSolve::solve: a.first > b.first)
      auto gt = [](const E& a, const E& b) { return a.first > b.first; };
      std::priority_queue<E, std::vector<E>, decltype(gt)> pq(gt);
      std::vector<E> popped;
      for (int i = 0; i < n; ++i) {
        pq.push(base[i]);
        if (i % 2 == 1) { popped.push_back(pq.top()); pq.pop(); }
      }
      while (!pq.empty()) { popped.push_back(pq.top()); pq.pop(); }
      out("pq", n, code, 0, popped);
    }
  }
  return 0;
}
