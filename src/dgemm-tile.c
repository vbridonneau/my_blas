#include <assert.h>
#include "dgemm.h"
#include "algonum.h"
#include <cblas.h>

#define BLOCK_SIZE 100
void my_dgemm_tile(const CBLAS_LAYOUT Order, 
      const enum CBLAS_TRANSPOSE TransA, 
      const enum CBLAS_TRANSPOSE TransB, 
      const int M, 
      const int N, 
      const int K, 
      const double alpha, 
      const double **A, 
      const int lda, 
      const double **B, 
      const int ldb, 
      const double beta, 
      double **C, 
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

	double * blockC = C[j * BLOCK_SIZE + i];

      for(int l = 0; l < i_blk_size; ++l) {
        for(int c = 0; c < j_blk_size; ++c) {
          blockC[l + c*BLOCK_SIZE] *= beta;
        }
      }

      for(int k = 0; k < kbloc; ++k) {
        int k_blk_size = (k < kbloc - 1 || !lastk) ? BLOCK_SIZE : lastk;
        if(!k_blk_size || !i_blk_size || !j_blk_size) continue;
        int shiftA     = (transA) * (k + i*lda) + (1-transA) * (i + k*lda);
        int shiftB     = (transB) * (j + k*ldb) + (1-transB) * (k + j*ldb);

	double * blockA = (transA)? A[k*BLOCK_SIZE + i] : A[i*BLOCK_SIZE + k];
	double * blockB = (transB)? B[j*BLOCK_SIZE + k] : A[k*BLOCK_SIZE + j];

        my_dgemm_scal_openmp(Order, TransA, TransB,
			  i_blk_size,
			  j_blk_size,
			  k_blk_size,
			  alpha,
			  blockA, BLOCK_SIZE,
			  blockB, BLOCK_SIZE,
			  1., blockC, ldc);
      }
    }
  }
}
