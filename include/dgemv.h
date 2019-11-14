#ifndef DEF_DGEMV_H
#define DEF_DGEMV_H
#include "definition.h"

void my_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
	      const int M, const int N, const double alpha, const double *A,
	      const int lda, const double *X, const int incX, const double beta,
	      double *Y, const int incY);

#endif//DEF_DGEMV_H
