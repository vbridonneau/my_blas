#ifndef DEF_DTRSM_H
#define DEF_DTRSM_H
#include "algonum.h"

void my_dtrsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE Side, const CBLAS_UPLO Uplo, CBLAS_TRANSPOSE transA, const CBLAS_DIAG Diag, const int M, const int N, const double alpha, const double *A, const int lda, double * B, const int ldb);
#endif//DEF_DGETRF_H
