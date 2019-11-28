#ifndef DEF_DGETRF_H
#define DEF_DGETRF_H

#include "algonum.h"

void my_dgetf2(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda );

void my_dgetrf(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda );

void my_dgetrf_omp(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda );

void my_pdgetrf(const CBLAS_LAYOUT Order, int M, int N, double* A, int lda );

#endif//DEF_DGETRF_H
