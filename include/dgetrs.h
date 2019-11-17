#ifndef DEF_DGETRS_H
#define DEF_DGETRS_H
#include "algonum.h"

int my_dgetrs(const CBLAS_TRANSPOSE trans, const int n, const int nrhs, const double * a, const int lda, const double * ipiv, double * b, const int ldb);

#endif//DEF_DGETRS_H
