#include "dger.h"

void my_dger(const int M, const int N, const double alpha,
    const double *X, const int incX, const double *Y,
    const int incY, double *A, const int lda) {
    for(int i = 0, stepX = 0; i < M; ++i, stepX += incX) {
        double tmp = X[stepX];
        for (int j = 0, stepY = 0; j < N; ++j, stepY += incY) {
            A[i + j*lda] += alpha*Y[stepY]*tmp;
        }
    }
}

