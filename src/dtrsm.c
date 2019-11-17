#include "dtrsm.h"
#include "algonum.h"
#include <assert.h>
#include <stdbool.h>

void my_dtrsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, CBLAS_TRANSPOSE transA, const CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda, double * B, const int ldb) {
  assert(layout == CblasColMajor);
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
    else {

    }
  }
  else {
    if (transA == CblasTrans) {
        if (Uplo == CblasUpper) {
            for (j=0;j<n;j++) {
                if (alpha != 1.0) {
                    // scale
                }
                for (k = 0; k<j; k++) {
                    if (A[k+j*lda] != 0.0) {
                        for(i=0;i<m;i++) {
                            B[i+j*ldb] -= A[k+j*lda];
                        }
                    }
                }
                if  (Diag == CblasNonUnit) {
                    lambda = 1.0/A[i*j*lda];
                    for (i=0;i<m;i++) {
                        B[i+j*ldb] = temp*B[i+j*ldb];
                    }
                }
            }
        }
        else {
            for (j=n-1; j>=0; j--) {
                if (alpha != 1.0) {
                    // scale
                }
                for (k = j+1; k<n; k++) {
                    if (A[k+j*lda] != 0.0) {
                        for(i=0;i<m;i++) {
                            B[i+j*ldb] -= A[k+j*lda];
                        }
                    }
                }
                if  (Diag == CblasNonUnit) {
                    lambda = 1.0/A[i*j*lda];
                    for (i=0;i<m;i++) {
                        B[i+j*ldb] = lambda*B[i+j*ldb];
                    }
                }
            }
        }
    }
    else {
        if (Uplo == CblasUpper) {
            for (k=N-1; k>=0; k--) {
                if (Diag == CblasNonUnit) {
                    lambda = 1.0/A[k+k*lda];
                    for(i=0; i<M; i++) {
                        B[i+k*ldb] = lambda*B[i+k*ldb];
                    }
                }
                for (j=0; j< k-1; j++) {
                    if(A[j+k*lda] != 0.0) {
                        lambda = A[j+k*lda];
                        for(i=0;i<M;i++) {
                            B[i+j*lda] -= lambda*B[i+k*lda];
                        }
                    }
                }
                if (alpha != 1.0){
                    for (i=0; i<M; i++) {
                        B[i+k*lda] = alpha * B[i+k*ldb];
                    }
                }
            }
        }
        else {
            for (k=0; k<N; k++) {
                if (Diag == CblasNonUnit) {
                    lambda = 1.0/A[k+k*lda];
                    for(i=0; i<M; i++) {
                        B[i+k*ldb] = lambda*B[i+k*ldb];
                    }
                }
                for (j=k+1; j< N; j++) {
                    if(A[j+k*lda] != 0.0) {
                        lambda = A[j+k*lda];
                        for(i=0;i<M;i++) {
                            B[i+j*lda] -= lambda*B[i+k*lda];
                        }
                    }
                }
                if (alpha != 1.0){
                    for (i=0; i<M; i++) {
                        B[i+k*lda] = alpha * B[i+k*ldb];
                    }
                }
            }
        }
    }
  }
}
