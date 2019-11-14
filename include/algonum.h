#ifndef _algonum_h_
#define _algonum_h_

#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include "flops.h"
#include "perf.h"

/**
 * Testing function for the matrix-matrix product
 */
typedef void (*dgemm_fct_t)( CBLAS_LAYOUT layout,
			     CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
			     int M, int N, int K,
			     double alpha, const double *A, int lda,
			                   const double *B, int ldb,
			     double beta,        double *C, int ldc );

int testone_dgemm( dgemm_fct_t dgemm,
                   CBLAS_TRANSPOSE transA,
                   CBLAS_TRANSPOSE transB,
                   int M, int N, int K, int check );
int testall_dgemm( dgemm_fct_t tested_dgemm );

/**
 * Testing function for the LU factorization
 */
typedef void (*dgetrf_fct_t)( CBLAS_LAYOUT layout,
	          	      int m, int n, double *a, int lda );

int testone_dgetrf( dgetrf_fct_t dgetrf,
                    int M, int N, int check );
int testall_dgetrf( dgetrf_fct_t tested_dgetrf );


/**
 * Function that you need to implement
 */
void my_dgemm_seq( CBLAS_LAYOUT layout,
                   CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                   int M, int N, int K,
                   double alpha, const double *A, int lda,
                                 const double *B, int ldb,
                   double beta,        double *C, int ldc );

void my_dgemm_scal_openmp( CBLAS_LAYOUT layout,
                           CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
                           int M, int N, int K,
                           double alpha, const double *A, int lda,
                                         const double *B, int ldb,
                           double beta,        double *C, int ldc );

void my_dgemm_bloc_openmp( CBLAS_LAYOUT layout,
                           CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
                           int M, int N, int K,
                           double alpha, const double *A, int lda,
                                         const double *B, int ldb,
                           double beta,        double *C, int ldc );

void my_dtrsm_openmp( CBLAS_LAYOUT layout,
                      CBLAS_SIDE side, CBLAS_UPLO uplo,
                      CBLAS_TRANSPOSE transA, CBLAS_DIAG diag,
                      int M, int N,
                      double alpha, const double *A, int lda,
	                                  double *B, int ldb );

void my_dgetrf_seq( CBLAS_LAYOUT layout,
                    int m, int n, double *a, int lda );

void my_dgetrf_openmp( CBLAS_LAYOUT layout,
                       int m, int n, double *a, int lda );

/**
 * Helpers for trhe conversion from lapack layout to tile layout
 */
double ** lapack2tile( int M, int N, int b, const double *Alapack, int lda );
void      tile2lapack( int M, int N, int b, const double **Atile, double *A, int lda );
void      tileFree( int M, int N, int b, double **A );

/**
 * Testing function for the tiled matrix-matrix product
 */
typedef void (*dgemm_tiled_fct_t)( CBLAS_LAYOUT layout,
                                   CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                                   int M, int N, int K, int b,
                                   double alpha, const double **A,
                                                 const double **B,
                                   double beta,        double **C );

int testone_dgemm_tiled( dgemm_tiled_fct_t dgemm,
                         CBLAS_TRANSPOSE transA,
                         CBLAS_TRANSPOSE transB,
                         int M, int N, int K, int b, int check );
int testall_dgemm_tiled( dgemm_tiled_fct_t tested_dgemm );

/**
 * Testing function for the LU factorization
 */
typedef void (*dgetrf_tiled_fct_t)( CBLAS_LAYOUT layout,
                                    int m, int n, int b, double **a );

int testone_dgetrf_tiled( dgetrf_tiled_fct_t dgetrf,
                          int M, int N, int b, int check );
int testall_dgetrf_tiled( dgetrf_tiled_fct_t tested_dgetrf );


/**
 * Function that you need to implement
 */
void my_dgemm_tiled_openmp( CBLAS_LAYOUT layout,
                            CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
                            int M, int N, int K, int b,
                            double alpha, const double **A,
                                          const double **B,
                            double beta,        double **C );

void my_dgetrf_tiled_openmp( CBLAS_LAYOUT layout,
                             int m, int n, int b, double **a );

#endif /* _algonum_h_ */
