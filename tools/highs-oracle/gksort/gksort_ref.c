/* Reference outputs for the GKlib sorts behind libmetis isortd / ikvsorti.
 *
 * Compiles the HiGHS-vendored GK_MKQSORT template (LGPL-2.1+, see
 * GKlib/gk_mksort.h) only into this test tool; nothing from it is copied into
 * the Rust crate. The Rust side is an independent implementation written from a
 * description of the sort (crates/rusolve-ekk/src/ipm/hipo/orderings/metis/gklib.rs),
 * and this tool's output is the evidence that the two order ties identically.
 *
 * Build and run (see gen.sh):
 *   cc -O1 -I$HIGHS/extern/metis/GKlib gksort_ref.c -o gksort_ref
 *   ./gksort_ref > crates/rusolve-ekk/tests/fixtures/gksort/gksort.txt
 *
 * Output, one case per line:
 *   d  <n> <in_0> ... <in_n-1> | <out_0> ... <out_n-1>      isortd
 *   kv <n> <key_0> ... <key_n-1> | <val_0> ... <val_n-1>    ikvsorti, val_i = i on input
 *   kvh <n> <seed> <keyrange> | <fnv1a64 of out (key,val) i32 LE>  large ikvsorti
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "gk_mksort.h"

typedef int32_t idx_t;
typedef struct {
  idx_t key;
  idx_t val;
} ikv_t;

static void isortd(size_t n, idx_t *base) {
#define i_gt(a, b) ((*a) > (*b))
  GK_MKQSORT(idx_t, base, n, i_gt);
#undef i_gt
}

static void ikvsorti(size_t n, ikv_t *base) {
#define ikey_lt(a, b) ((a)->key < (b)->key)
  GK_MKQSORT(ikv_t, base, n, ikey_lt);
#undef ikey_lt
}

/* xorshift32; the Rust test does not regenerate inputs, it reads them. */
static uint32_t rng = 2463534242u;
static uint32_t next(void) {
  rng ^= rng << 13;
  rng ^= rng >> 17;
  rng ^= rng << 5;
  return rng;
}

static idx_t draw(uint32_t range) { return (idx_t)(next() % range) - (idx_t)(range / 3); }

static void case_d(size_t n, uint32_t range) {
  idx_t *a = malloc((n ? n : 1) * sizeof *a);
  for (size_t i = 0; i < n; i++) a[i] = draw(range);
  printf("d %zu", n);
  for (size_t i = 0; i < n; i++) printf(" %d", a[i]);
  isortd(n, a);
  printf(" |");
  for (size_t i = 0; i < n; i++) printf(" %d", a[i]);
  printf("\n");
  free(a);
}

static void fill_kv(ikv_t *a, size_t n, uint32_t range, int shape) {
  for (size_t i = 0; i < n; i++) {
    idx_t k;
    switch (shape) {
      case 1: k = (idx_t)i / 3; break;                    /* ascending runs */
      case 2: k = (idx_t)(n - i) / 2; break;              /* descending runs */
      case 3: k = (i % 7 == 0) ? draw(range) : 5; break;  /* mostly equal */
      case 4: k = (idx_t)(i % 4); break;                  /* sawtooth */
      default: k = draw(range); break;
    }
    a[i].key = k;
    a[i].val = (idx_t)i;
  }
}

static void case_kv(size_t n, uint32_t range, int shape) {
  ikv_t *a = malloc((n ? n : 1) * sizeof *a);
  fill_kv(a, n, range, shape);
  printf("kv %zu", n);
  for (size_t i = 0; i < n; i++) printf(" %d", a[i].key);
  ikvsorti(n, a);
  printf(" |");
  for (size_t i = 0; i < n; i++) printf(" %d", a[i].val);
  printf("\n");
  free(a);
}

/* Large case: inputs come from a second generator the Rust test reimplements. */
static uint32_t lcg;
static idx_t lkey(uint32_t range) {
  lcg = lcg * 1103515245u + 12345u;
  return (idx_t)((lcg >> 8) % range);
}

static void case_kvh(size_t n, uint32_t seed, uint32_t range) {
  ikv_t *a = malloc(n * sizeof *a);
  lcg = seed;
  for (size_t i = 0; i < n; i++) {
    a[i].key = lkey(range);
    a[i].val = (idx_t)i;
  }
  ikvsorti(n, a);
  uint64_t h = 1469598103934665603ull;
  for (size_t i = 0; i < n; i++) {
    int32_t w[2] = {a[i].key, a[i].val};
    const unsigned char *p = (const unsigned char *)w;
    for (int b = 0; b < 8; b++) {
      h ^= p[b];
      h *= 1099511628211ull;
    }
  }
  printf("kvh %zu %u %u | %llu\n", n, seed, range, (unsigned long long)h);
  free(a);
}

int main(void) {
  static const size_t sizes[] = {0, 1, 2, 3, 7, 8, 9, 10, 11, 16, 17, 18, 19, 25, 31, 40, 64, 97, 150, 300};
  static const uint32_t ranges[] = {1, 2, 3, 5, 17, 1000};
  const size_t ns = sizeof sizes / sizeof *sizes, nr = sizeof ranges / sizeof *ranges;
  for (size_t s = 0; s < ns; s++)
    for (size_t r = 0; r < nr; r++) case_d(sizes[s], ranges[r]);
  for (size_t s = 0; s < ns; s++) {
    for (size_t r = 0; r < nr; r++) case_kv(sizes[s], ranges[r], 0);
    for (int shape = 1; shape <= 4; shape++) case_kv(sizes[s], 9, shape);
  }
  static const uint32_t big_ranges[] = {2, 10, 1000, 1000000};
  for (uint32_t seed = 1; seed <= 3; seed++)
    for (size_t r = 0; r < 4; r++) case_kvh(20000, seed, big_ranges[r]);
  return 0;
}
