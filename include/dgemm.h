#ifndef DGEMM_DEF_H
#define DGEMM_DEF_H
#include "algonum.h"

/* Sequential versions */
void my_dgemm_scalaire(const CBLAS_LAYOUT Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);

void my_dgemm(const CBLAS_LAYOUT Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);

/* OpenMP versions */
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
			  const int ldc);

void my_dgemm_omp(const CBLAS_LAYOUT Order, 
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
			  const int ldc);

/* Tiled versions */
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
		   double **C);

#endif//DGEMM_DEF_H
