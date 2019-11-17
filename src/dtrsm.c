#include "dtrsm.h"
#include "algonum.h"
#include <assert.h>
#include <stdbool.h>

void my_dtrsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, CBLAS_TRANSPOSE transA, const CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda, double * B, const int ldb) {
  assert(layout == CblasColMajor);
  /* assert(side == 'l'); */
  /* assert(uplo == 'u'); */
  double lambda;

  if (M == 0 || N == 0) return;

  /* scale 0. */
  if (alpha == 0.) {
    for(int j = 0; j < N; ++j) {
      for (int i = 0; i < M; ++i) {
	B[i + j * ldb] = 0.;
      }
    }
    return;
  }

  // FIXME: a ne devrait pas etre modifié
  // TODO: prendre en compte les parametres side et diag

  /* Left side : X * op( A ) = alpha * B */
  if (Side == CblasLeft) {
    /* B = alpha * inv(A ** t) * B */
    if (transA == CblasTrans) {
      /* A is a lower triangular */
      if (Uplo == CblasLower) {
	for(int j = 0; j < N; ++j) {
	  for(int i = M - 1; i >= 0; --i) {
	    lambda = alpha * B[i + j*ldb];
	    for(int k = i + 1; k < M; ++k) {
	      lambda -= B[k + j*ldb] * A[k + i*lda];
	    }
	    /* The diagonal is A[i + i*lda] (Otherwise : 1.) */
	    /* Relevent when solving A = L*U as we use A to store
	     both L and U, so Diag(L) is full of 1. . */
	    if (Diag == CblasNonUnit) lambda /= A[i*(1 + lda)];
	    B[i + j*ldb] = lambda;
	  }
	}
      }
      /* A is triangular upper */
      else if (Uplo == CblasUpper) {
	for(int j= 0; j < N; ++j) {
	  for(int i = 0; i < M; ++i) {
	    lambda = alpha * B[i + j*ldb];
	    for(int k = 0; k < i; ++k) {
	      lambda -=  A[k + i*lda] * B[k + j*ldb];
	    }
	    if (Diag == CblasNonUnit) lambda /= A[i*(1 + lda)];
	    B[i + j*ldb] = lambda;
	  }
	}
      }
    }
  }
  else {
    assert(false);
  }
}
