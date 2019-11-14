#include "dgemv.h"
#include "dscal.h"
#include "ddot.h"
#include <assert.h>
#include <string.h>

static void dgemv_transA_ColMajor(const int M, const int N, const double alpha,
				  const double *A, const int lda,
				  const double *X, const int incX,
				  const double beta, double *Y,
				  const int incY) {
  for(int i = 0, stepY = 0; i < N; ++i, stepY += incY) {
    double dotAX = alpha*my_ddot(M, A + i*lda, 1, X, 1);
    Y[stepY] += dotAX;
  }
}

static void dgemv_A_ColMajor(const int M, const int N, const double alpha,
			     const double *A, const int lda,
			     const double *X, const int incX,
			     const double beta, double *Y,
			     const int incY) {
  for(int j = 0, stepX = 0; j < N; ++j, stepX += incX) {
    for(int i = 0, stepY = 0; i < M; ++i, stepY += incY) {
      Y[stepY] += alpha * A[i + j*lda] * X[stepX];;
    }
  }
}

void my_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
	      const int M, const int N, const double alpha, const double *A,
	      const int lda, const double *X, const int incX, const double beta,
	      double *Y, const int incY) {
  assert(Order == COLUMN_MAJOR);
  assert(Trans == 'T' || Trans == 't' || Trans == 'N' || Trans == 'n');
  int transA = (Trans == 't' || Trans == 'T');
  int lenY = (transA) ? N : M;
  /* Scale Y before all */
  if (beta != 1.) {
    if (beta == 0. && incY == 1) {
      memset (Y, 0, sizeof(double) * lenY);
    } else {
      my_dscal(lenY, beta, Y, incY);
    }
  }

  if (transA) {
    dgemv_transA_ColMajor(M, N, alpha, A, lda, X, incX, 1., Y, incY);
  } else {
    dgemv_A_ColMajor(M, N, alpha, A, lda, X, incX, 1., Y, incY);
  }
}
