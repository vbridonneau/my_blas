#ifndef DEF_DGESV_H
#define DEF_DGESV_H

#include "algonum.h"

void dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int *ldb, int *info);

#endif//DEF_DGETRF_H
