#include <assert.h>
#include "dgemm.h"
#include "algonum.h"
#include <cblas.h>
#include <stdio.h>

#define BLOCK_SIZE 100
void my_dgemm_tile(const CBLAS_LAYOUT Order, 
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

  for (int i = 0; i < MT; ++i) {
    int i_blk_size = (i < MT - 1 || !lastm) ? b : lastm;
    for (int j = 0; j < NT; ++j) {
      int j_blk_size = (j < NT - 1 || !lastn) ? b : lastn;

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

        my_dgemm_scalaire(Order, TransA, TransB,
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
