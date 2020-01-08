#include "algonum.h"
#include <assert.h>

void my_dgemm_seq( CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
                   CBLAS_TRANSPOSE transB, const int M, const int N,
                   const int K, const double alpha, const double *A,
                   const int lda, const double *B, const int ldb,
                   const double beta, double *C, const int ldc )
{
  int m, n, k;

  if ( transA == CblasNoTrans ) {
    if ( transB == CblasNoTrans ) {
      for( k=0; k<K; k++ ) {
	double lbeta = k == 0 ? beta : 1.;
	for( m=0; m<M; m++ ) {
	  for( n=0; n<N; n++ ) {
	    C[ ldc * n + m ] = alpha * A[ lda * k + m ] * B[ ldb * n + k ]
	      +              beta * C[ ldc * n + m ];
	  }
	}
      }
    }
    else {
      for( k=0; k<K; k++ ) {
	double lbeta = k == 0 ? beta : 1.;
	for( m=0; m<M; m++ ) {
	  for( n=0; n<N; n++ ) {
	    C[ ldc * n + m ] = alpha * A[ lda * k + m ] * B[ ldb * k + n ]
	      +              beta * C[ ldc * n + m ];
	  }
	}
      }
    }
  }
  else {
    if ( transB == CblasNoTrans ) {
      for( k=0; k<K; k++ ) {
	double lbeta = k == 0 ? beta : 1.;
	for( m=0; m<M; m++ ) {
	  for( n=0; n<N; n++ ) {
	    C[ ldc * n + m ] = alpha * A[ lda * m + k ] * B[ ldb * n + k ]
	      +              beta * C[ ldc * n + m ];
	  }
	}
      }
    }
    else {
      for( k=0; k<K; k++ ) {
	double lbeta = k == 0 ? beta : 1.;
	for( m=0; m<M; m++ ) {
	  for( n=0; n<N; n++ ) {
	    C[ ldc * n + m ] = alpha * A[ lda * m + k ] * B[ ldb * k + n ]
	      +              beta * C[ ldc * n + m ];
	  }
	}
      }
    }
  }
}

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
  if (layout != CblasColMajor) return;
  if ((TransA == CblasTrans) && (TransB == CblasNoTrans)){
    my_dgemm_tAB(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  } else if (TransA == CblasNoTrans && TransB == CblasNoTrans) {
    my_dgemm_AB(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  } else if (TransA == CblasNoTrans && TransB == CblasTrans) {
    my_dgemm_AtB(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  } else if (TransA == CblasTrans && TransB == CblasTrans) {
    my_dgemm_tAtB(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
  } else {
    assert(0);
  }
}



#define BLOCK_SIZE 100
void my_dgemm_bloc_openmp(const CBLAS_LAYOUT Order, 
      const enum CBLAS_TRANSPOSE TransA, 
      const enum CBLAS_TRANSPOSE TransB, 
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
    #pragma omp parallel for
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

        my_dgemm_scal_openmp(Order, TransA, TransB,
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
void my_dgemm_tiled_openmp(const CBLAS_LAYOUT Order, 
		   const enum CBLAS_TRANSPOSE TransA, 
		   const enum CBLAS_TRANSPOSE TransB, 
		   const int M, 
		   const int N, 
		   const int K,
		   const int b,
		   const double alpha, 
		   const double **A, 
		   const double **B, 
		   const double beta, 
		   double **C) {
  if (Order != CblasColMajor) return;
  int transA = (TransA == CblasTrans);
  int transB = (TransB == CblasTrans);

  if(b == 0) return;

  /* Upper bound for tile in m, n, k */
  int MT, NT, KT;
  MT = (M + b - 1) / b;
  NT = (N + b - 1) / b;
  KT = (K + b - 1) / b;

  
  /* Reminders */
  int lastk = K % b;
  int lastm = M % b;
  int lastn = N % b;

  for (int j = 0; j < NT; ++j) {
    int j_blk_size = (j < NT - 1 || !lastn) ? b : lastn;

    #pragma omp parallel for
    for (int i = 0; i < MT; ++i) {
      int i_blk_size = (i < MT - 1 || !lastm) ? b : lastm;


      double * blockC = C[j * MT + i];

      for(int l = 0; l < i_blk_size; ++l) {
        for(int c = 0; c < j_blk_size; ++c) {
          blockC[l + c*b] *= beta;
        }
      }

      for(int k = 0; k < KT; ++k) {
        int k_blk_size = (k < KT - 1 || !lastk) ? b : lastk;
        if(!k_blk_size || !i_blk_size || !j_blk_size) continue;

	const double * blockA = (!transA)? A[k*MT + i] : A[i*KT + k];
	const double * blockB = (!transB)? B[j*KT + k] : B[k*NT + j];

        my_dgemm_seq(Order, TransA, TransB,
			  i_blk_size,
			  j_blk_size,
			  k_blk_size,
			  alpha,
			  blockA, b,
			  blockB, b,
			  1., blockC, b);
      }
    }
  }
}

/* To make sure we use the right prototype */
static dgemm_fct_t valig_mygemm __attribute__ ((unused)) = my_dgemm_seq;
