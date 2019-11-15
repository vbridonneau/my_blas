#ifndef DEF_DGETRF_H
#define DEF_DGETRF_H

#include "algonum.h"

int my_dgetf2(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda );

int my_dgetrf(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda );

#endif//DEF_DGETRF_H
