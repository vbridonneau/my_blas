#ifndef DEF_DGER_H
#define DEF_DGER_H

void my_dger(const int M, const int N, const double alpha,
    const double *X, const int incX, const double *Y,
    const int incY, double *A, const int lda);

#endif//DEF_DGER_H
