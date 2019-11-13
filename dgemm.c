#include "dgemm.h"
#include <math.h>
#include <assert.h>
#include <ctype.h> // tolower && toupper

/* Transposed A */
static inline void my_dgemm_tAB_kij(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int k = 0; k < K; ++k) { /* row of B, col of tA */
    for (int i = 0; i < M; ++i) { /* row of tA */
      for (int j = 0; j < N; ++j) { /* col of B */
	C[i + j*ldc] += alpha * A[k + i*lda] * B[k + j*ldb];
      }
    }
  }
}

static inline void my_dgemm_tAB_ijk(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int i = 0; i < M; ++i) { /* row of tA */
    for (int j = 0; j < N; ++j) { /* col of B */
      for (int k = 0; k < K; ++k) { /* row of B, col of tA */
	C[i + j*ldc] += alpha * A[k + i*lda] * B[k + j*ldb];
      }
    }
  }
}

static inline void my_dgemm_tAB_jik(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
    for (int i = 0; i < M; ++i) { /* row of tA */
      for (int k = 0; k < K; ++k) { /* row of B, col of tA */
	C[i + j*ldc] += alpha * A[k + i*lda] * B[k + j*ldb];
      }
    }
  }
}

/* Non transposed A */
static inline void my_dgemm_AB_kij(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int k = 0; k < K; ++k) { /* row of B, col of A */
    for (int i = 0; i < M; ++i) { /* row of A */
      for (int j = 0; j < N; ++j) { /* col of B */
	C[i + j*ldc] += alpha * A[i + k*lda] * B[k + j*ldb];
      }
    }
  }
}


static inline void my_dgemm_AB_jik(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int j = 0; j < N; ++j) { /* col of B */
    for (int i = 0; i < M; ++i) { /* row of A */
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	C[i + j*ldc] += alpha * A[i + k*lda] * B[k + j*ldb];
      }
    }
  }
}

static inline void my_dgemm_AB_ijk(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  for (int i = 0; i < M; ++i) { /* row of A */
    for (int j = 0; j < N; ++j) { /* col of B */
      for (int k = 0; k < K; ++k) { /* row of B, col of A */
	C[i + j*ldc] += alpha * A[i + k*lda] * B[k + j*ldb];
      }
    }
  }
}

/* Scalar Version */
void my_dgemm_scalaire(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
    if (Order != COLUMN_MAJOR) return;
    if ((toupper(TransA) == 'T') && (toupper(TransB) == 'N')){
      my_dgemm_tAB_kij(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    } else if (toupper(TransA) == 'N' && toupper(TransB) == 'N') {
      my_dgemm_AB_kij(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    }
}

#define BLOCK_SIZE 128
void my_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc) {
  if (Order != COLUMN_MAJOR) return;
  if ((toupper(TransA) != 'N' && toupper(TransA) != 'T')
      || toupper(TransB) != 'N') return;
  int transA = (toupper(TransA) == 'T');
  int kbloc  = (int)ceil((double)K/(double)BLOCK_SIZE), lastk = K % BLOCK_SIZE;
  int mbloc  = (int)ceil((double)M/(double)BLOCK_SIZE), lastm = M % BLOCK_SIZE;
  int nbloc  = (int)ceil((double)N/(double)BLOCK_SIZE), lastn = N % BLOCK_SIZE;
  for(int k = 0; k < kbloc; ++k) {
    for (int i = 0; i < mbloc; ++i) {
      for (int j = 0; j < nbloc; ++j) {
	int k_blk_size = (k < kbloc - 1 || !lastk) ? BLOCK_SIZE : lastk;
	int i_blk_size = (i < mbloc - 1 || !lastm) ? BLOCK_SIZE : lastm;
	int j_blk_size = (j < nbloc - 1 || !lastn) ? BLOCK_SIZE : lastn;
	if(!k_blk_size || !i_blk_size || !j_blk_size) continue;
	int shiftA     = (transA) ? k + i*lda: i + k*lda;
	my_dgemm_scalaire(Order, TransA, TransB, i_blk_size, j_blk_size, k_blk_size, alpha, A + shiftA*BLOCK_SIZE, lda, B + (k + j*ldb)*BLOCK_SIZE, ldb, beta, C + (i + j*ldc)*BLOCK_SIZE, ldc);
      }
    }
  }
}
