#ifndef DEF_DGETRF_H
#define DEF_DGETRF_H

int my_dgetf2(const enum CBLAS_ORDER Order, int m, int n, double* a, int lda, int* ipiv );

int my_dgetrf(const enum CBLAS_ORDER Order, int m, int n, double* a, int lda, int* ipiv );

#endif//DEF_DGETRF_H
