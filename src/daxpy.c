#include "daxpy.h"

void my_daxpy(const int N, const double a, const double *X, const int incX,
	   double *Y, const int incY) {
  int stepX = 0, stepY = 0;
  for(int i = 0; i < N; ++ i, stepX += incX, stepY += incY) {
    Y[stepY] += a * X[stepX];
  }
}

