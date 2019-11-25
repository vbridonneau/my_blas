#include "dgemm.h"
#include "algonum.h"

static inline void my_dgemm_tAB(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
#pragma omp parallel for
    for (int i = 0; i < M; ++i) { /* row of tA */
      C[i + j*ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of tA */
	C[i + j*ldc] += alpha * A[k + i*lda] * B[k + j*ldb];
      }
    }
  }
}

static inline void my_dgemm_AB(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
#pragma omp parallel for
    for (int i = 0; i < M; ++i) { /* row of A */
      C[i + j * ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	C[i + j*ldc] += alpha * A[i + k*lda] * B[k + j*ldb];
      }
    }
  }
}

static inline void my_dgemm_AtB(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for(int j = 0; j < N; ++j) {
#pragma omp parallel for
    for(int i = 0; i < M; ++i) {
      C[i + j*ldc] *= beta;
    }
  }

  for (int k = 0; k < K; ++k) { /* row of B, col of A */
    for (int j = 0; j < N; ++j) { /* col of B */
#pragma omp parallel for
      for (int i = 0; i < M; ++i) { /* row of A */
	C[i + j*ldc] += alpha * A[i + k*lda] * B[j + k*ldb];
	  /* This avoid conditionnal jump, but add computing complexity */
	  // + (bool_k * 1. + (!bool_k) * beta)*C[i + j*ldc];
      }
    }
  }
}

static inline void my_dgemm_tAtB(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
#pragma omp parallel for
    for (int i = 0; i < M; ++i) { /* row of A */
      C[i + j * ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	C[i + j*ldc] += alpha * A[k + i*lda] * B[j + k*ldb];
      }
    }
  }
}

void my_dgemm_scal_openmp(const CBLAS_LAYOUT layout,
                          const CBLAS_TRANSPOSE TransA,
			  const CBLAS_TRANSPOSE TransB,
			  const int M,
			  const int N,
			  const int K,
			  const double alpha,
			  const double *A,
			  const int lda,
			  const double *B,
			  const int ldb,
			  const double beta,
			  double *C,
			  const int ldc) {
  if (Order != CblasColMajor) return;
  if ((TransA == CblasTrans) && (TransB == CblasNoTrans)){
    my_dgemm_tAB_ijk(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  } else if (TransA == CblasNoTrans && TransB == CblasNoTrans) {
    my_dgemm_AB_jik(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  } else if (TransA == CblasNoTrans && TransB == CblasTrans) {
    my_dgemm_AtB_kji(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  } else if (TransA == CblasTrans && TransB == CblasTrans) {
    my_dgemm_tAtB_ijk(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  } else {
    assert(0);
  }
}
