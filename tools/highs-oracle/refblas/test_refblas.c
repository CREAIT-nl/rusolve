/* refblas unit checks against hand-computed cases (exact: every input and
 * result is a small dyadic rational, so == is the right comparison).
 * MIT License, Copyright (c) 2026 the rusolve authors.
 *
 * Covers every routine and every (uplo, trans, diag, side) case HiPO uses,
 * plus strided and negative increments. Exit 0 iff all pass. */
#include <stdio.h>

#include "refblas.h"

static int failures = 0;

static void check(const char* name, const double* got, const double* want,
                  int n) {
  int i;
  for (i = 0; i < n; i++) {
    if (got[i] != want[i]) {
      fprintf(stderr, "FAIL %s [%d]: got %.17g want %.17g\n", name, i, got[i],
              want[i]);
      failures++;
      return;
    }
  }
  printf("ok   %s\n", name);
}

int main(void) {
  /* ----- level 1 */
  {
    double x[3] = {1, 2, 3}, y[3] = {1, 1, 1}, w[3] = {3, 5, 7};
    cblas_daxpy(3, 2.0, x, 1, y, 1);
    check("daxpy unit", y, w, 3);
  }
  {
    double x[3] = {1, 2, 3}, z[4] = {0, 0, 0, 0}, w[4] = {1, 0, 0, 3};
    cblas_daxpy(2, 1.0, x, 2, z, 3);
    check("daxpy strided", z, w, 4);
  }
  {
    /* incx=-1: x is visited from its last element */
    double x[3] = {1, 2, 3}, y[3] = {0, 0, 0}, w[3] = {3, 2, 1};
    cblas_daxpy(3, 1.0, x, -1, y, 1);
    check("daxpy negative inc", y, w, 3);
  }
  {
    double x[2] = {1, 2}, y[2] = {5, 6}, w[2] = {5, 6};
    cblas_daxpy(2, 0.0, x, 1, y, 1);
    check("daxpy alpha=0 no-op", y, w, 2);
  }
  {
    double x[3] = {1, 2, 3}, y[6] = {0, 0, 0, 0, 0, 0},
           w[6] = {1, 0, 2, 0, 3, 0};
    cblas_dcopy(3, x, 1, y, 2);
    check("dcopy strided", y, w, 6);
  }
  {
    double x[4] = {1, 2, 3, 4}, w[4] = {0.5, 2, 1.5, 4};
    cblas_dscal(2, 0.5, x, 2);
    check("dscal strided", x, w, 4);
  }
  {
    double x[2] = {1, 2}, w[2] = {1, 2};
    cblas_dscal(2, 3.0, x, -1); /* netlib: incx<=0 returns */
    check("dscal negative inc no-op", x, w, 2);
  }
  {
    double x[3] = {1, 2, 3}, y[3] = {4, 5, 6}, wx[3] = {4, 5, 6},
           wy[3] = {1, 2, 3};
    cblas_dswap(3, x, 1, y, 1);
    check("dswap x", x, wx, 3);
    check("dswap y", y, wy, 3);
  }
  {
    /* swap columns 0 and 2 of a 2x3 col-major array (lda=2), as HiPO's
     * swapCols does on one buffer */
    double a[6] = {1, 2, 3, 4, 5, 6}, w[6] = {5, 6, 3, 4, 1, 2};
    cblas_dswap(2, a, 1, a + 4, 1);
    check("dswap within one array", a, w, 6);
  }
  {
    /* swap rows 0 and 1 of a 2x2 col-major array: stride lda=2 */
    double a[4] = {1, 2, 3, 4}, w[4] = {2, 1, 4, 3};
    cblas_dswap(2, a, 2, a + 1, 2);
    check("dswap strided", a, w, 4);
  }

  /* ----- level 2. A = [1 3; 2 4] col-major */
  {
    double a[4] = {1, 2, 3, 4}, x[2] = {1, 1}, y[2] = {10, 10},
           w[2] = {4 * 1 + 10 * 0.5, 6 * 1 + 10 * 0.5};
    /* y := 1*A*x + 0.5*y = [4,6] + [5,5] */
    cblas_dgemv(CblasColMajor, CblasNoTrans, 2, 2, 1.0, a, 2, x, 1, 0.5, y, 1);
    check("dgemv N", y, w, 2);
  }
  {
    double a[4] = {1, 2, 3, 4}, x[2] = {1, 1}, y[2] = {10, 10},
           w[2] = {-3, -7};
    /* y := -1*A^T*x + 0*y = -[3,7] */
    cblas_dgemv(CblasColMajor, CblasTrans, 2, 2, -1.0, a, 2, x, 1, 0.0, y, 1);
    check("dgemv T beta=0", y, w, 2);
  }
  {
    /* 3x2 A = [1 4; 2 5; 3 6], lda=3, x strided by 2, y by -1 */
    double a[6] = {1, 2, 3, 4, 5, 6}, x[3] = {1, 99, 2}, y[3] = {0, 0, 0},
           w[3] = {15, 12, 9};
    /* A*x = [1+8, 2+10, 3+12] = [9,12,15]; incy=-1 stores element i at
     * y[2-i] */
    cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 2, 1.0, a, 3, x, 2, 0.0, y,
                -1);
    check("dgemv N strided/negative", y, w, 3);
  }
  {
    /* U = [2 1; 0 4], solve U x = b with b = U*[1,1] = [3,4] */
    double a[4] = {2, 0, 1, 4}, x[2] = {3, 4}, w[2] = {1, 1};
    cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 2, a, 2,
                x, 1);
    check("dtrsv U N nonunit", x, w, 2);
  }
  {
    /* unit upper U = [1 3; 0 1]: U x = [7,2] -> x = [1,2] */
    double a[4] = {9, 0, 3, 9}, x[2] = {7, 2}, w[2] = {1, 2};
    cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, 2, a, 2, x,
                1);
    check("dtrsv U N unit", x, w, 2);
  }
  {
    /* unit upper U^T = [1 0; 3 1]: U^T x = [1,5] -> x = [1,2] */
    double a[4] = {9, 0, 3, 9}, x[2] = {1, 5}, w[2] = {1, 2};
    cblas_dtrsv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, 2, a, 2, x,
                1);
    check("dtrsv U T unit", x, w, 2);
  }
  {
    /* L = [1 0; 2 1] unit: L x = [1,4] -> [1,2] */
    double a[4] = {1, 2, 0, 1}, x[2] = {1, 4}, w[2] = {1, 2};
    cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasUnit, 2, a, 2, x,
                1);
    check("dtrsv L N unit", x, w, 2);
  }
  {
    /* L = [2 0; 1 4] nonunit, L^T = [2 1; 0 4]: L^T x = [3,4] -> [1,1];
     * x strided by 2 */
    double a[4] = {2, 1, 0, 4}, x[3] = {3, 99, 4}, w[3] = {1, 99, 1};
    cblas_dtrsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, 2, a, 2,
                x, 2);
    check("dtrsv L T nonunit strided", x, w, 3);
  }
  {
    /* packed upper U = [2 1; 0 4] -> ap = {2, 1, 4} */
    double ap[3] = {2, 1, 4}, x[2] = {3, 4}, w[2] = {1, 1};
    cblas_dtpsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 2, ap,
                x, 1);
    check("dtpsv U N", x, w, 2);
  }
  {
    /* U^T = [2 0; 1 4]: U^T x = [2,5] -> [1,1] */
    double ap[3] = {2, 1, 4}, x[2] = {2, 5}, w[2] = {1, 1};
    cblas_dtpsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, 2, ap, x,
                1);
    check("dtpsv U T", x, w, 2);
  }
  {
    /* packed lower L = [2 0 0; 1 4 0; 1 1 2] -> ap = {2,1,1, 4,1, 2};
     * L*[1,1,1] = [2,5,4] */
    double ap[6] = {2, 1, 1, 4, 1, 2}, x[3] = {2, 5, 4}, w[3] = {1, 1, 1};
    cblas_dtpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, 3, ap,
                x, 1);
    check("dtpsv L N", x, w, 3);
  }
  {
    /* L^T = [2 1 1; 0 4 1; 0 0 2]; L^T*[1,1,1] = [4,5,2] */
    double ap[6] = {2, 1, 1, 4, 1, 2}, x[3] = {4, 5, 2}, w[3] = {1, 1, 1};
    cblas_dtpsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, 3, ap, x,
                1);
    check("dtpsv L T", x, w, 3);
  }
  {
    /* A (2x2 zero) += 2 * x y^T, x = [1,2], y = [3,0] (y(2)=0 skipped) */
    double a[4] = {1, 1, 1, 1}, x[2] = {1, 2}, y[2] = {3, 0},
           w[4] = {7, 13, 1, 1};
    cblas_dger(CblasColMajor, 2, 2, 2.0, x, 1, y, 1, a, 2);
    check("dger", a, w, 4);
  }

  /* ----- level 3 */
  {
    /* C = A * B^T with B = I */
    double a[4] = {1, 2, 3, 4}, b[4] = {1, 0, 0, 1}, c[4] = {0, 0, 0, 0},
           w[4] = {1, 2, 3, 4};
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 2, 2, 2, 1.0, a, 2, b,
                2, 0.0, c, 2);
    check("dgemm N T", c, w, 4);
  }
  {
    /* C := -1*A^T*B + 1*C, A = [1 3; 2 4], B = [1 0; 1 1], C = ones.
     * A^T = [1 2; 3 4]; A^T B = [3 2; 7 4] */
    double a[4] = {1, 2, 3, 4}, b[4] = {1, 1, 0, 1}, c[4] = {1, 1, 1, 1},
           w[4] = {-2, -6, -1, -3};
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 2, 2, 2, -1.0, a, 2, b,
                2, 1.0, c, 2);
    check("dgemm T N beta=1", c, w, 4);
  }
  {
    /* C := A*B + 2*C, A 2x3 = [1 2 3; 4 5 6], B 3x1 = [1,1,1], C = [1,1] */
    double a[6] = {1, 4, 2, 5, 3, 6}, b[3] = {1, 1, 1}, c[2] = {1, 1},
           w[2] = {8, 17};
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 2, 1, 3, 1.0, a, 2,
                b, 3, 2.0, c, 2);
    check("dgemm N N", c, w, 2);
  }
  {
    /* C := A^T*B^T, A = [1 3; 2 4], B = I: C = A^T = [1 2; 3 4] */
    double a[4] = {1, 2, 3, 4}, b[4] = {1, 0, 0, 1}, c[4] = {9, 9, 9, 9},
           w[4] = {1, 3, 2, 4};
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, 2, 2, 2, 1.0, a, 2, b, 2,
                0.0, c, 2);
    check("dgemm T T", c, w, 4);
  }
  {
    /* lower C := A A^T, A = [1 2; 3 4]: A A^T = [5 11; 11 25]; the upper
     * entry keeps its 9 */
    double a[4] = {1, 3, 2, 4}, c[4] = {0, 0, 9, 0}, w[4] = {5, 11, 9, 25};
    cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, 2, 2, 1.0, a, 2, 0.0,
                c, 2);
    check("dsyrk L N", c, w, 4);
  }
  {
    /* upper C := A^T A + C, A = [1 2; 3 4]: A^T A = [10 14; 14 20] */
    double a[4] = {1, 3, 2, 4}, c[4] = {1, 9, 1, 1}, w[4] = {11, 9, 15, 21};
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, 2, 2, 1.0, a, 2, 1.0, c,
                2);
    check("dsyrk U T", c, w, 4);
  }
  {
    /* HiPO's case: left, upper, trans, unit. U = [1 2; 0 1],
     * U^T = [1 0; 2 1]; solve U^T X = B, B = U^T*[1 1; 1 2] = [1 1; 3 4] */
    double a[4] = {7, 0, 2, 7}, b[4] = {1, 3, 1, 4}, w[4] = {1, 1, 1, 2};
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasUnit, 2,
                2, 1.0, a, 2, b, 2);
    check("dtrsm L U T unit", b, w, 4);
  }
  {
    /* left upper notrans nonunit with alpha=2: U = [2 1; 0 4],
     * X = U^-1 (2B), B = [3,4]/2*... use B = [1.5, 2] -> 2B = [3,4] -> [1,1] */
    double a[4] = {2, 0, 1, 4}, b[2] = {1.5, 2}, w[2] = {1, 1};
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
                CblasNonUnit, 2, 1, 2.0, a, 2, b, 2);
    check("dtrsm L U N nonunit alpha", b, w, 2);
  }
  {
    /* left lower notrans: L = [2 0; 1 4], L*[1,1] = [2,5] */
    double a[4] = {2, 1, 0, 4}, b[2] = {2, 5}, w[2] = {1, 1};
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                CblasNonUnit, 2, 1, 1.0, a, 2, b, 2);
    check("dtrsm L L N", b, w, 2);
  }
  {
    /* left lower trans: L^T = [2 1; 0 4], L^T*[1,1] = [3,4] */
    double a[4] = {2, 1, 0, 4}, b[2] = {3, 4}, w[2] = {1, 1};
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit,
                2, 1, 1.0, a, 2, b, 2);
    check("dtrsm L L T", b, w, 2);
  }
  {
    /* right upper notrans: X U = B, U = [2 1; 0 4], [1 1] U = [2 5] */
    double a[4] = {2, 0, 1, 4}, b[2] = {2, 5}, w[2] = {1, 1};
    cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans,
                CblasNonUnit, 1, 2, 1.0, a, 2, b, 1);
    check("dtrsm R U N", b, w, 2);
  }
  {
    /* right lower notrans: X L = B, L = [2 0; 1 4], [1 1] L = [3 4] */
    double a[4] = {2, 1, 0, 4}, b[2] = {3, 4}, w[2] = {1, 1};
    cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasNoTrans,
                CblasNonUnit, 1, 2, 1.0, a, 2, b, 1);
    check("dtrsm R L N", b, w, 2);
  }
  {
    /* right upper trans: X U^T = B, U^T = [2 0; 1 4], [1 1] U^T = [3 4] */
    double a[4] = {2, 0, 1, 4}, b[2] = {3, 4}, w[2] = {1, 1};
    cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasTrans,
                CblasNonUnit, 1, 2, 1.0, a, 2, b, 1);
    check("dtrsm R U T", b, w, 2);
  }
  {
    /* right lower trans: X L^T = B, L^T = [2 1; 0 4], [1 1] L^T = [2 5] */
    double a[4] = {2, 1, 0, 4}, b[2] = {2, 5}, w[2] = {1, 1};
    cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans,
                CblasNonUnit, 1, 2, 1.0, a, 2, b, 1);
    check("dtrsm R L T", b, w, 2);
  }

  if (failures) {
    fprintf(stderr, "refblas: %d check(s) failed\n", failures);
    return 1;
  }
  printf("refblas: all checks passed\n");
  return 0;
}
