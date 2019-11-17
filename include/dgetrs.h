#ifndef DEF_DGETRS_H
#define DEF_DGETRS_H

void dgetrs(
    char trans,
    int n,
    int nrhs,
    double * a,
    int lda,
    double * ipiv,
    double * b,
    int ldb,
    int info
);

#endif//DEF_DGETRS_H
