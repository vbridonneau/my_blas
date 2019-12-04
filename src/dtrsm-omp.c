#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include "cblas.h"
#include "dtrsm.h"
#include "algonum.h"

/*
layout  CHARACTER*1
              = CblasColMajor   storage order is column
              = CblasRowMajor   storage order is row
SIDE    CHARACTER*1
              = CblasLeft     op( A )*X = alpha*B.
              = CblasRight    X*op( A ) = alpha*B.
UPLO    CHARACTER*1
              = CblasUpper  A is an upper triangular matrix.
              = CblasLower   A is a lower triangular matrix.
TRANSA  CHARACTER*1
              = CblasNoTrans => op( A ) = A
              = CblasTrans => op( A ) = A**T
DIAG    CHARACTER*1
              = CblasUnit if A is unit triangular
              = CblasNonUnit if A is not unit triangular
M       number of rows of B
N       INTEGER number of columns of B
ALPHA   scalar
A       DOUBLE PRECISION array LDA*k
           where k is m when SIDE = 'L'
             and k is n when SIDE = 'R'
LDA     dimension of A
B       DOUBLE PRECISION array LDB*N
LDB	    dimension of B
*/
void my_dtrsm_omp(const CBLAS_LAYOUT layout, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, CBLAS_TRANSPOSE transA, const CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda, double * B, const int ldb) {
  assert(layout == CblasColMajor);
  double lambda;

  if (M == 0 || N == 0) return;

  /* scale 0. */
  if (alpha == 0.) {
    for(int j = 0; j < N; ++j) {
      memset(B + j * ldb, 0, M * sizeof(double));
    }
    return;
  }

  /* Left side : op( A ) * X = alpha * B */
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
	    /* The diagonal is A[i + i*lda] (Otherwise : 1.) */
	    if (Diag == CblasNonUnit) lambda /= A[i*(1 + lda)];
	    B[i + j*ldb] = lambda;
	  }
	}
      }
    }
    /* B = alpha * inv(A) * B */
    else {
      /* A is triangular Upper */
      if (Uplo == CblasUpper) {
	for (int j = 0; j < N; ++j) {
	  if (alpha != 1.) {
	    #pragma omp parallel for static(1)
	    for (int i = 0; i < M; i++) {
	      B[i + j*ldb] *= alpha;
	    }
	  }
	  for (int k = M - 1; k >= 0; --k) {
	    if (B[k + j*ldb]) {
	      if (Diag == CblasNonUnit) B[k + j*ldb] /= A[k*(1 + lda)];
	      lambda = B[k + j*ldb];
	      for (int i = 0; i < k; ++i) {
		B[i + j*ldb] -= lambda*A[i + k*lda];
	      }
	    }
	  }
	}
      }
      /* A is lower triangular */
      else {
	for (int j = 0; j < N; ++j) {
	  #pragma omp parallel for static(1)
	  for (int i = 0; i < M; i++) {
	    B[i + j*ldb] *= alpha;
	  }
	  for (int k = 0; k < M; ++k) {
	    if (B[k + j*ldb] != 0.) {
	      if (Diag == CblasNonUnit)
		B[k + j*ldb] /= A[k*(1 + lda)];
	      lambda = B[k + j*ldb];
	      for (int i = k+1; i < M; ++i) {
		B[i + j*ldb] -= lambda*A[i + k*lda];
	      }
	    }
	  }
	}
      }
    }
  }
  /* Right side : X * op( A ) = alpha*B */
  else {
    /* X = alpha * B * inv(A) */
    if (transA == CblasNoTrans) {
      /* A is upper triangular */
      if (Uplo == CblasUpper) {
	for (int j = 0; j < N;j++) {
	  if (alpha != 1.0) {
	    #pragma omp parallel for static(1)
	    for (int i = 0; i < M; ++i) {
	      B[i + j*ldb] *= alpha;
	    }
	  }
	  for (int k = 0; k < j - 1; k++) {
	    if (A[k+j*lda] != 0.0) {
	      for(int i = 0;i < M;i++) {
		B[i + j*ldb] -= A[k+j*lda] * B[i + k*ldb];
	      }
	    }
	  }
	  if  (Diag == CblasNonUnit) {
	    lambda = 1.0/A[j*(1 + lda)];
	    #pragma omp parallel for static(1)
	    for (int i = 0;i < M;i++) {
	      B[i+j*ldb] = lambda*B[i+j*ldb];
	    }
	  }
	}
      }
      /* A is lower triangular */
      else {
	for (int j = N-1; j>=0; --j) {
	  if (alpha != 1.0) {
	    #pragma omp parallel for static(1)
	    for (int i = 0; i < M; ++i) {
	      B[i + j*ldb] *= alpha;
	    }
	  }
	  for (int k = j+1; k < N; ++k) {
	    if (A[k+j*lda] != 0.0) {
	      for(int i = 0;i < M; ++i) {
		B[i + j*ldb] -= A[k + j*lda] * B[i + k*ldb];
	      }
	    }
	  }
	  if  (Diag == CblasNonUnit) {
	    lambda = 1.0/A[j*(1 + lda)];
	    #pragma omp parallel for static(1)
	    for (int i = 0;i < M;i++) {
	      B[i+j*ldb] = lambda*B[i+j*ldb];
	    }
	  }
	}
      }
    }
    /* X = alpha * B * inv(A ** t) */
    else {
      /* A is upper triangular */
      if (Uplo == CblasUpper) {
	for (int k = N-1; k >= 0; --k) {
	  if (Diag == CblasNonUnit) {
	    lambda = 1.0/A[k+k*lda];
	    #pragma omp parallel for static(1)
	    for(int i = 0; i < M; i++) {
	      B[i+k*ldb] = lambda*B[i+k*ldb];
	    }
	  }
	  for (int j = 0; j < k; ++j) {
	    if(A[j+k*lda] != 0.0) {
	      lambda = A[j+k*lda];
	      for(int i = 0;i < M; ++i) {
		B[i+j*ldb] -= lambda*B[i+k*ldb];
	      }
	    }
	  }
	  if (alpha != 1.0){
	    #pragma omp parallel for static(1)
	    for (int i = 0; i < M; i++) {
	      B[i+k*ldb] = alpha * B[i+k*ldb];
	    }
	  }
	}
      }
      /* A is lower triangular */
      else {
	for (int k = 0; k < N; ++k) {
	  if (Diag == CblasNonUnit) {
	    lambda = 1.0/A[k+k*lda];
	    #pragma omp parallel for static(1)
	    for(int i = 0; i < M; ++i) {
	      B[i+k*ldb] = lambda*B[i+k*ldb];
	    }
	  }
	  for (int j = k+1; j < N; j++) {
	    if(A[j+k*lda] != 0.0) {
	      lambda = A[j+k*lda];
	      for(int i = 0;i < M;i++) {
		B[i+j*ldb] -= lambda*B[i+k*ldb];
	      }
	    }
	  }
	  if (alpha != 1.0) {
	    #pragma omp parallel for static(1)
	    for (int i = 0; i < M; ++i) {
	      B[i+k*lda] = alpha * B[i+k*ldb];
	    }
	  }
	}
      }
    }
  }
}
