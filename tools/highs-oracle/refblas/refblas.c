/* refblas: a reference CBLAS for the rusolve HiPO oracle.
 * MIT License, Copyright (c) 2026 the rusolve authors.
 *
 * Why it exists: HiPO's results depend on the BLAS it links, and an
 * optimised BLAS (OpenBLAS, Accelerate) blocks, unrolls and fuses its loops
 * differently per CPU. The oracle links this file instead, so every fixture
 * is produced by the loop order written here, and the Rust port mirrors it
 * statement by statement (crates/rusolve-ekk/src/ipm/hipo/blas.rs).
 *
 * Loop order: netlib reference BLAS (LAPACK 3.12 BLAS/SRC/d*.f), translated
 * to 0-based C. Deliberate simplifications that do not change any floating
 * point result:
 *  - the netlib level-1 unrolled unit-stride loops (daxpy by 4, dcopy by 7,
 *    dscal by 5, dswap by 3) are written as the plain loop: every element is
 *    computed by the same single expression either way;
 *  - netlib's separate unit-stride and general-stride branches of the level
 *    2 routines are merged into the general-stride form, which visits the
 *    same elements with the same operations in the same order;
 *  - argument checking (XERBLA) is not done: HiPO only passes valid
 *    arguments.
 * Everything else (the early returns, the beta==0 / beta!=1 tests, the
 * x(j)!=0 / y(j)!=0 / b(k,j)!=0 skips, the temp accumulators and the order of
 * operands in each expression) is netlib's.
 *
 * Column-major only. CblasConjTrans is CblasTrans (real data). A
 * CblasRowMajor call aborts: HiPO never makes one, and a silent wrong answer
 * would be worse than a crash.
 *
 * Build: -O2 -ffp-contract=off (no FMA contraction), C99. */
#include "refblas.h"

#include <stdio.h>
#include <stdlib.h>

#define A_(i, j) a[(size_t)(i) + (size_t)(j) * (size_t)lda]
#define B_(i, j) b[(size_t)(i) + (size_t)(j) * (size_t)ldb]
#define C_(i, j) c[(size_t)(i) + (size_t)(j) * (size_t)ldc]

static void refblas_col_major_only(const char* name, int order) {
  if (order != CblasColMajor) {
    fprintf(stderr, "refblas: %s called with order %d (only ColMajor)\n",
            name, order);
    abort();
  }
}

/* start index of a strided vector of length n (netlib KX, 0-based) */
static blasint refblas_start(blasint n, blasint inc) {
  return inc > 0 ? 0 : -(n - 1) * inc;
}

/* ---------------------------------------------------------------- level 1 */

void cblas_daxpy(const blasint n, const double alpha, const double* x,
                 const blasint incx, double* y, const blasint incy) {
  blasint i, ix, iy;
  if (n <= 0) return;
  if (alpha == 0.0) return;
  ix = refblas_start(n, incx);
  iy = refblas_start(n, incy);
  for (i = 0; i < n; i++) {
    y[iy] = y[iy] + alpha * x[ix];
    ix += incx;
    iy += incy;
  }
}

void cblas_dcopy(const blasint n, const double* x, const blasint incx,
                 double* y, const blasint incy) {
  blasint i, ix, iy;
  if (n <= 0) return;
  ix = refblas_start(n, incx);
  iy = refblas_start(n, incy);
  for (i = 0; i < n; i++) {
    y[iy] = x[ix];
    ix += incx;
    iy += incy;
  }
}

void cblas_dscal(const blasint n, const double alpha, double* x,
                 const blasint incx) {
  blasint i, ix;
  if (n <= 0 || incx <= 0 || alpha == 1.0) return;
  ix = 0;
  for (i = 0; i < n; i++) {
    x[ix] = alpha * x[ix];
    ix += incx;
  }
}

void cblas_dswap(const blasint n, double* x, const blasint incx, double* y,
                 const blasint incy) {
  blasint i, ix, iy;
  double temp;
  if (n <= 0) return;
  ix = refblas_start(n, incx);
  iy = refblas_start(n, incy);
  for (i = 0; i < n; i++) {
    temp = x[ix];
    x[ix] = y[iy];
    y[iy] = temp;
    ix += incx;
    iy += incy;
  }
}

/* ---------------------------------------------------------------- level 2 */

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa, const blasint M,
                 const blasint n, const double alpha, const double* a,
                 const blasint lda, const double* x, const blasint incx,
                 const double beta, double* y, const blasint incy) {
  const blasint m = M;
  blasint i, j, lenx, leny, kx, ky, ix, iy, jx, jy;
  double temp;
  refblas_col_major_only("cblas_dgemv", order);
  if (m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0)) return;
  if (transa == CblasNoTrans) {
    lenx = n;
    leny = m;
  } else {
    lenx = m;
    leny = n;
  }
  kx = refblas_start(lenx, incx);
  ky = refblas_start(leny, incy);
  /* y := beta*y */
  if (beta != 1.0) {
    iy = ky;
    if (beta == 0.0) {
      for (i = 0; i < leny; i++) {
        y[iy] = 0.0;
        iy += incy;
      }
    } else {
      for (i = 0; i < leny; i++) {
        y[iy] = beta * y[iy];
        iy += incy;
      }
    }
  }
  if (alpha == 0.0) return;
  if (transa == CblasNoTrans) {
    /* y := alpha*A*x + y */
    jx = kx;
    for (j = 0; j < n; j++) {
      temp = alpha * x[jx];
      iy = ky;
      for (i = 0; i < m; i++) {
        y[iy] = y[iy] + temp * A_(i, j);
        iy += incy;
      }
      jx += incx;
    }
  } else {
    /* y := alpha*A**T*x + y */
    jy = ky;
    for (j = 0; j < n; j++) {
      temp = 0.0;
      ix = kx;
      for (i = 0; i < m; i++) {
        temp = temp + A_(i, j) * x[ix];
        ix += incx;
      }
      y[jy] = y[jy] + alpha * temp;
      jy += incy;
    }
  }
}

void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const blasint n, const double* ap, double* x,
                 const blasint incx) {
  blasint i, j, k, kk, kx, ix, jx;
  double temp;
  const int nounit = diag == CblasNonUnit;
  refblas_col_major_only("cblas_dtpsv", order);
  if (n == 0) return;
  kx = refblas_start(n, incx);
  /* kk is the 0-based packed index of the current diagonal element */
  if (transa == CblasNoTrans) {
    /* x := inv(A)*x */
    if (uplo == CblasUpper) {
      kk = (n * (n + 1)) / 2 - 1;
      jx = kx + (n - 1) * incx;
      for (j = n - 1; j >= 0; j--) {
        if (x[jx] != 0.0) {
          if (nounit) x[jx] = x[jx] / ap[kk];
          temp = x[jx];
          ix = jx;
          for (k = kk - 1; k >= kk - j; k--) {
            ix -= incx;
            x[ix] = x[ix] - temp * ap[k];
          }
        }
        jx -= incx;
        kk -= j + 1;
      }
    } else {
      kk = 0;
      jx = kx;
      for (j = 0; j < n; j++) {
        if (x[jx] != 0.0) {
          if (nounit) x[jx] = x[jx] / ap[kk];
          temp = x[jx];
          ix = jx;
          for (k = kk + 1; k <= kk + n - 1 - j; k++) {
            ix += incx;
            x[ix] = x[ix] - temp * ap[k];
          }
        }
        jx += incx;
        kk += n - j;
      }
    }
  } else {
    /* x := inv(A**T)*x */
    if (uplo == CblasUpper) {
      kk = 0; /* packed index of the first element of column j */
      jx = kx;
      for (j = 0; j < n; j++) {
        temp = x[jx];
        ix = kx;
        for (k = kk; k <= kk + j - 1; k++) {
          temp = temp - ap[k] * x[ix];
          ix += incx;
        }
        if (nounit) temp = temp / ap[kk + j];
        x[jx] = temp;
        jx += incx;
        kk += j + 1;
      }
    } else {
      kk = (n * (n + 1)) / 2 - 1; /* packed index of the last element of col j */
      kx = kx + (n - 1) * incx;
      jx = kx;
      for (j = n - 1; j >= 0; j--) {
        temp = x[jx];
        ix = kx;
        for (k = kk; k >= kk - (n - 1 - j) + 1; k--) {
          temp = temp - ap[k] * x[ix];
          ix -= incx;
        }
        if (nounit) temp = temp / ap[kk - (n - 1 - j)];
        x[jx] = temp;
        jx -= incx;
        kk -= n - j;
      }
    }
  }
  (void)i;
}

void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const blasint n, const double* a, const blasint lda, double* x,
                 const blasint incx) {
  blasint i, j, kx, ix, jx;
  double temp;
  const int nounit = diag == CblasNonUnit;
  refblas_col_major_only("cblas_dtrsv", order);
  if (n == 0) return;
  kx = refblas_start(n, incx);
  if (transa == CblasNoTrans) {
    /* x := inv(A)*x */
    if (uplo == CblasUpper) {
      jx = kx + (n - 1) * incx;
      for (j = n - 1; j >= 0; j--) {
        if (x[jx] != 0.0) {
          if (nounit) x[jx] = x[jx] / A_(j, j);
          temp = x[jx];
          ix = jx;
          for (i = j - 1; i >= 0; i--) {
            ix -= incx;
            x[ix] = x[ix] - temp * A_(i, j);
          }
        }
        jx -= incx;
      }
    } else {
      jx = kx;
      for (j = 0; j < n; j++) {
        if (x[jx] != 0.0) {
          if (nounit) x[jx] = x[jx] / A_(j, j);
          temp = x[jx];
          ix = jx;
          for (i = j + 1; i < n; i++) {
            ix += incx;
            x[ix] = x[ix] - temp * A_(i, j);
          }
        }
        jx += incx;
      }
    }
  } else {
    /* x := inv(A**T)*x */
    if (uplo == CblasUpper) {
      jx = kx;
      for (j = 0; j < n; j++) {
        temp = x[jx];
        ix = kx;
        for (i = 0; i < j; i++) {
          temp = temp - A_(i, j) * x[ix];
          ix += incx;
        }
        if (nounit) temp = temp / A_(j, j);
        x[jx] = temp;
        jx += incx;
      }
    } else {
      kx = kx + (n - 1) * incx;
      jx = kx;
      for (j = n - 1; j >= 0; j--) {
        temp = x[jx];
        ix = kx;
        for (i = n - 1; i > j; i--) {
          temp = temp - A_(i, j) * x[ix];
          ix -= incx;
        }
        if (nounit) temp = temp / A_(j, j);
        x[jx] = temp;
        jx -= incx;
      }
    }
  }
}

void cblas_dger(const enum CBLAS_ORDER order, const blasint m, const blasint n,
                const double alpha, const double* x, const blasint incx,
                const double* y, const blasint incy, double* a,
                const blasint lda) {
  blasint i, j, kx, ix, jy;
  double temp;
  refblas_col_major_only("cblas_dger", order);
  if (m == 0 || n == 0 || alpha == 0.0) return;
  jy = refblas_start(n, incy);
  kx = refblas_start(m, incx);
  for (j = 0; j < n; j++) {
    if (y[jy] != 0.0) {
      temp = alpha * y[jy];
      ix = kx;
      for (i = 0; i < m; i++) {
        A_(i, j) = A_(i, j) + x[ix] * temp;
        ix += incx;
      }
    }
    jy += incy;
  }
}

/* ---------------------------------------------------------------- level 3 */

void cblas_dgemm(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_TRANSPOSE transb, const blasint m,
                 const blasint n, const blasint k, const double alpha,
                 const double* a, const blasint lda, const double* b,
                 const blasint ldb, const double beta, double* c,
                 const blasint ldc) {
  blasint i, j, l;
  double temp;
  const int nota = transa == CblasNoTrans;
  const int notb = transb == CblasNoTrans;
  refblas_col_major_only("cblas_dgemm", order);
  if (m == 0 || n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0)) return;
  if (alpha == 0.0) {
    if (beta == 0.0) {
      for (j = 0; j < n; j++)
        for (i = 0; i < m; i++) C_(i, j) = 0.0;
    } else {
      for (j = 0; j < n; j++)
        for (i = 0; i < m; i++) C_(i, j) = beta * C_(i, j);
    }
    return;
  }
  if (notb) {
    if (nota) {
      /* C := alpha*A*B + beta*C */
      for (j = 0; j < n; j++) {
        if (beta == 0.0) {
          for (i = 0; i < m; i++) C_(i, j) = 0.0;
        } else if (beta != 1.0) {
          for (i = 0; i < m; i++) C_(i, j) = beta * C_(i, j);
        }
        for (l = 0; l < k; l++) {
          temp = alpha * B_(l, j);
          for (i = 0; i < m; i++) C_(i, j) = C_(i, j) + temp * A_(i, l);
        }
      }
    } else {
      /* C := alpha*A**T*B + beta*C */
      for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
          temp = 0.0;
          for (l = 0; l < k; l++) temp = temp + A_(l, i) * B_(l, j);
          if (beta == 0.0)
            C_(i, j) = alpha * temp;
          else
            C_(i, j) = alpha * temp + beta * C_(i, j);
        }
      }
    }
  } else {
    if (nota) {
      /* C := alpha*A*B**T + beta*C */
      for (j = 0; j < n; j++) {
        if (beta == 0.0) {
          for (i = 0; i < m; i++) C_(i, j) = 0.0;
        } else if (beta != 1.0) {
          for (i = 0; i < m; i++) C_(i, j) = beta * C_(i, j);
        }
        for (l = 0; l < k; l++) {
          temp = alpha * B_(j, l);
          for (i = 0; i < m; i++) C_(i, j) = C_(i, j) + temp * A_(i, l);
        }
      }
    } else {
      /* C := alpha*A**T*B**T + beta*C */
      for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
          temp = 0.0;
          for (l = 0; l < k; l++) temp = temp + A_(l, i) * B_(j, l);
          if (beta == 0.0)
            C_(i, j) = alpha * temp;
          else
            C_(i, j) = alpha * temp + beta * C_(i, j);
        }
      }
    }
  }
}

void cblas_dsyrk(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans, const blasint n,
                 const blasint k, const double alpha, const double* a,
                 const blasint lda, const double beta, double* c,
                 const blasint ldc) {
  blasint i, j, l, i0, i1;
  double temp;
  const int upper = uplo == CblasUpper;
  refblas_col_major_only("cblas_dsyrk", order);
  if (n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0)) return;
  /* the stored triangle of column j is rows [i0, i1) */
#define REFBLAS_TRI(j_) \
  do {                  \
    if (upper) {        \
      i0 = 0;           \
      i1 = (j_) + 1;    \
    } else {            \
      i0 = (j_);        \
      i1 = n;           \
    }                   \
  } while (0)
  if (alpha == 0.0) {
    for (j = 0; j < n; j++) {
      REFBLAS_TRI(j);
      if (beta == 0.0) {
        for (i = i0; i < i1; i++) C_(i, j) = 0.0;
      } else {
        for (i = i0; i < i1; i++) C_(i, j) = beta * C_(i, j);
      }
    }
    return;
  }
  if (trans == CblasNoTrans) {
    /* C := alpha*A*A**T + beta*C */
    for (j = 0; j < n; j++) {
      REFBLAS_TRI(j);
      if (beta == 0.0) {
        for (i = i0; i < i1; i++) C_(i, j) = 0.0;
      } else if (beta != 1.0) {
        for (i = i0; i < i1; i++) C_(i, j) = beta * C_(i, j);
      }
      for (l = 0; l < k; l++) {
        if (A_(j, l) != 0.0) {
          temp = alpha * A_(j, l);
          for (i = i0; i < i1; i++) C_(i, j) = C_(i, j) + temp * A_(i, l);
        }
      }
    }
  } else {
    /* C := alpha*A**T*A + beta*C */
    for (j = 0; j < n; j++) {
      REFBLAS_TRI(j);
      for (i = i0; i < i1; i++) {
        temp = 0.0;
        for (l = 0; l < k; l++) temp = temp + A_(l, i) * A_(l, j);
        if (beta == 0.0)
          C_(i, j) = alpha * temp;
        else
          C_(i, j) = alpha * temp + beta * C_(i, j);
      }
    }
  }
#undef REFBLAS_TRI
}

void cblas_dtrsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side,
                 const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_DIAG diag, const blasint m, const blasint n,
                 const double alpha, const double* a, const blasint lda,
                 double* b, const blasint ldb) {
  blasint i, j, k;
  double temp;
  const int lside = side == CblasLeft;
  const int upper = uplo == CblasUpper;
  const int nounit = diag == CblasNonUnit;
  refblas_col_major_only("cblas_dtrsm", order);
  if (m == 0 || n == 0) return;
  if (alpha == 0.0) {
    for (j = 0; j < n; j++)
      for (i = 0; i < m; i++) B_(i, j) = 0.0;
    return;
  }
  if (lside) {
    if (transa == CblasNoTrans) {
      /* B := alpha*inv(A)*B */
      if (upper) {
        for (j = 0; j < n; j++) {
          if (alpha != 1.0)
            for (i = 0; i < m; i++) B_(i, j) = alpha * B_(i, j);
          for (k = m - 1; k >= 0; k--) {
            if (B_(k, j) != 0.0) {
              if (nounit) B_(k, j) = B_(k, j) / A_(k, k);
              for (i = 0; i < k; i++) B_(i, j) = B_(i, j) - B_(k, j) * A_(i, k);
            }
          }
        }
      } else {
        for (j = 0; j < n; j++) {
          if (alpha != 1.0)
            for (i = 0; i < m; i++) B_(i, j) = alpha * B_(i, j);
          for (k = 0; k < m; k++) {
            if (B_(k, j) != 0.0) {
              if (nounit) B_(k, j) = B_(k, j) / A_(k, k);
              for (i = k + 1; i < m; i++)
                B_(i, j) = B_(i, j) - B_(k, j) * A_(i, k);
            }
          }
        }
      }
    } else {
      /* B := alpha*inv(A**T)*B */
      if (upper) {
        for (j = 0; j < n; j++) {
          for (i = 0; i < m; i++) {
            temp = alpha * B_(i, j);
            for (k = 0; k < i; k++) temp = temp - A_(k, i) * B_(k, j);
            if (nounit) temp = temp / A_(i, i);
            B_(i, j) = temp;
          }
        }
      } else {
        for (j = 0; j < n; j++) {
          for (i = m - 1; i >= 0; i--) {
            temp = alpha * B_(i, j);
            for (k = i + 1; k < m; k++) temp = temp - A_(k, i) * B_(k, j);
            if (nounit) temp = temp / A_(i, i);
            B_(i, j) = temp;
          }
        }
      }
    }
  } else {
    if (transa == CblasNoTrans) {
      /* B := alpha*B*inv(A) */
      if (upper) {
        for (j = 0; j < n; j++) {
          if (alpha != 1.0)
            for (i = 0; i < m; i++) B_(i, j) = alpha * B_(i, j);
          for (k = 0; k < j; k++) {
            if (A_(k, j) != 0.0) {
              for (i = 0; i < m; i++) B_(i, j) = B_(i, j) - A_(k, j) * B_(i, k);
            }
          }
          if (nounit) {
            temp = 1.0 / A_(j, j);
            for (i = 0; i < m; i++) B_(i, j) = temp * B_(i, j);
          }
        }
      } else {
        for (j = n - 1; j >= 0; j--) {
          if (alpha != 1.0)
            for (i = 0; i < m; i++) B_(i, j) = alpha * B_(i, j);
          for (k = j + 1; k < n; k++) {
            if (A_(k, j) != 0.0) {
              for (i = 0; i < m; i++) B_(i, j) = B_(i, j) - A_(k, j) * B_(i, k);
            }
          }
          if (nounit) {
            temp = 1.0 / A_(j, j);
            for (i = 0; i < m; i++) B_(i, j) = temp * B_(i, j);
          }
        }
      }
    } else {
      /* B := alpha*B*inv(A**T) */
      if (upper) {
        for (k = n - 1; k >= 0; k--) {
          if (nounit) {
            temp = 1.0 / A_(k, k);
            for (i = 0; i < m; i++) B_(i, k) = temp * B_(i, k);
          }
          for (j = 0; j < k; j++) {
            if (A_(j, k) != 0.0) {
              temp = A_(j, k);
              for (i = 0; i < m; i++) B_(i, j) = B_(i, j) - temp * B_(i, k);
            }
          }
          if (alpha != 1.0)
            for (i = 0; i < m; i++) B_(i, k) = alpha * B_(i, k);
        }
      } else {
        for (k = 0; k < n; k++) {
          if (nounit) {
            temp = 1.0 / A_(k, k);
            for (i = 0; i < m; i++) B_(i, k) = temp * B_(i, k);
          }
          for (j = k + 1; j < n; j++) {
            if (A_(j, k) != 0.0) {
              temp = A_(j, k);
              for (i = 0; i < m; i++) B_(i, j) = B_(i, j) - temp * B_(i, k);
            }
          }
          if (alpha != 1.0)
            for (i = 0; i < m; i++) B_(i, k) = alpha * B_(i, k);
        }
      }
    }
  }
}
