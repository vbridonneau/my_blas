#include "dgemm.h"
#include <math.h>
#include <assert.h>
#include <ctype.h> // tolower && toupper
#include "algonum.h"

/* Transposed A */
static inline void my_dgemm_tAB_kij(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int k = 0; k < K; ++k) { /* row of B, col of tA */
    int bool_k = !!k;
    for (int i = 0; i < M; ++i) { /* row of tA */
      for (int j = 0; j < N; ++j) { /* col of B */
	C[i + j*ldc] = alpha * A[k + i*lda] * B[k + j*ldb]
	  /* This avoid conditionnal jump, but add computing complexity */
	  + (bool_k * 1. + (!bool_k) * beta)*C[i + j*ldc];
      }
    }
  }
}

static inline void my_dgemm_tAB_ijk(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int i = 0; i < M; ++i) { /* row of tA */
    for (int j = 0; j < N; ++j) { /* col of B */
      C[i + j*ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of tA */
	C[i + j*ldc] += alpha * A[k + i*lda] * B[k + j*ldb];
      }
    }
  }
}

static inline void my_dgemm_tAB_jik(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
    for (int i = 0; i < M; ++i) { /* row of tA */
      C[i + j*ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of tA */
	C[i + j*ldc] += alpha * A[k + i*lda] * B[k + j*ldb];
      }
    }
  }
}

/* Non transposed A */
static inline void my_dgemm_AB_kij(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int k = 0; k < K; ++k) { /* row of B, col of A */
    int bool_k = !!k;
    for (int i = 0; i < M; ++i) { /* row of A */
      for (int j = 0; j < N; ++j) { /* col of B */
	C[i + j*ldc] = alpha * A[i + k*lda] * B[k + j*ldb]
	  /* This avoid conditionnal jump, but add computing complexity */
	  + (bool_k * 1. + (!bool_k) * beta)*C[i + j*ldc];
      }
    }
  }
}


static inline void my_dgemm_AB_jik(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
    for (int i = 0; i < M; ++i) { /* row of A */
      C[i + j * ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	C[i + j*ldc] += alpha * A[i + k*lda] * B[k + j*ldb];
      }
    }
  }
}

static inline void my_dgemm_AB_ijk(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int i = 0; i < M; ++i) { /* row of A */
    for (int j = 0; j < N; ++j) { /* col of B */
      C[i + j * ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
      	C[i + j*ldc] += alpha * A[i + k*lda] * B[k + j*ldb];
      }
    }
  }
}

/* No transposed A, Transposed B */
static inline void my_dgemm_AtB_kji(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int k = 0; k < K; ++k) { /* row of B, col of A */
    int bool_k = !!k;
    for (int j = 0; j < N; ++j) { /* col of B */
      for (int i = 0; i < M; ++i) { /* row of A */
	C[i + j*ldc] = alpha * A[i + k*lda] * B[j + k*ldb]
	  /* This avoid conditionnal jump, but add computing complexity */
	  + (bool_k * 1. + (!bool_k) * beta)*C[i + j*ldc];
      }
    }
  }
}

static inline void my_dgemm_AtB_jik(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
    for (int i = 0; i < M; ++i) { /* row of A */
      C[i + j * ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	C[i + j*ldc] += alpha * A[i + k*lda] * B[j + k*ldb];
      }
    }
  }
}

static inline void my_dgemm_AtB_ijk(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int i = 0; i < M; ++i) { /* row of A */
    for (int j = 0; j < N; ++j) { /* col of B */
      C[i + j * ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	C[i + j*ldc] += alpha * A[i + k*lda] * B[j + k*ldb];
      }
    }
  }
}

/* Transposed A, Transposed B */
static inline void my_dgemm_tAtB_kij(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int k = 0; k < K; ++k) { /* row of tB, col of tA */
    int bool_k = !!k;
    for (int j = 0; j < N; ++j) { /* col of tB */
      for (int i = 0; i < M; ++i) { /* row of tA */
	C[i + j*ldc] = alpha * A[k + i*lda] * B[j + k*ldb]
	  /* This avoid conditionnal jump, but add computing complexity */
	  + (bool_k * 1. + (!bool_k) * beta)*C[i + j*ldc];
      }
    }
  }
}

static inline void my_dgemm_tAtB_jik(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
    for (int i = 0; i < M; ++i) { /* row of A */
      C[i + j * ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	      C[i + j*ldc] += alpha * A[k + i*lda] * B[j + k*ldb];
      }
    }
  }
}

static inline void my_dgemm_tAtB_ijk(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int i = 0; i < M; ++i) { /* row of A */
    for (int j = 0; j < N; ++j) { /* col of B */
      C[i + j * ldc] *= beta;
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	      C[i + j*ldc] += alpha * A[k + i*lda] * B[j + k*ldb];
      }
    }
  }
}

/* Scalar Version */
void my_dgemm_scalaire(const CBLAS_LAYOUT Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
    if (Order != CblasColMajor) return;
    if ((TransA == CblasTrans) && (TransB == CblasNoTrans)){
      my_dgemm_tAB_ijk(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    } else if (TransA == CblasNoTrans && TransB == CblasNoTrans) {
      my_dgemm_AB_jik(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    } else if (TransA == CblasNoTrans && TransB == CblasTrans) {
      my_dgemm_AtB_jik(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    } else if (TransA == CblasTrans && TransB == CblasTrans) {
      my_dgemm_tAtB_ijk(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    } else {
      assert(0);
    }
}

#define BLOCK_SIZE 100
void my_dgemm(const CBLAS_LAYOUT Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  if (Order != CblasColMajor) return;
  int transA = (TransA == CblasTrans);
  int transB = (TransB == CblasTrans);
  int kbloc, nbloc, mbloc;

  if (K % BLOCK_SIZE)
    kbloc = (K + BLOCK_SIZE)/BLOCK_SIZE;
  else
    kbloc =  K/ BLOCK_SIZE;
  if (M % BLOCK_SIZE)
    mbloc = (M + BLOCK_SIZE)/BLOCK_SIZE;
  else
    mbloc = M / BLOCK_SIZE;
  if (N % BLOCK_SIZE)
    nbloc = (N + BLOCK_SIZE)/BLOCK_SIZE;
  else
    nbloc = N / BLOCK_SIZE;
  
  /* Reminders */
  int lastk = K % BLOCK_SIZE;
  int lastm = M % BLOCK_SIZE;
  int lastn = N % BLOCK_SIZE;

  for (int i = 0; i < mbloc; ++i) {
    int i_blk_size = (i < mbloc - 1 || !lastm) ? BLOCK_SIZE : lastm;
    for (int j = 0; j < nbloc; ++j) {
      int j_blk_size = (j < nbloc - 1 || !lastn) ? BLOCK_SIZE : lastn;

      double *Ctranslated = C + (i + j*ldc)*BLOCK_SIZE;
      for(int l = 0; l < i_blk_size; ++l) {
	for(int c = 0; c < j_blk_size; ++c) {
	  Ctranslated[l + c*ldc] *= beta;
	}
      }

      for(int k = 0; k < kbloc; ++k) {
	int k_blk_size = (k < kbloc - 1 || !lastk) ? BLOCK_SIZE : lastk;
	if(!k_blk_size || !i_blk_size || !j_blk_size) continue;
	int shiftA     = (transA) * (k + i*lda) + (1-transA) * (i + k*lda);
	int shiftB     = (transB) * (j + k*ldb) + (1-transB) * (k + j*ldb);
        my_dgemm_scalaire(Order, TransA, TransB,
			  i_blk_size,
			  j_blk_size,
			  k_blk_size,
			  alpha,
			  A + shiftA*BLOCK_SIZE, lda,
			  B + shiftB*BLOCK_SIZE, ldb,
			  1., Ctranslated, ldc);
      }
    }
  }
}
