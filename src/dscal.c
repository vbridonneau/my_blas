#include "dscal.h"

void my_dscal(const int N, const double alpha, double *X, const int incX) {
  for(int i = 0, stepX = 0; i < N; ++i, stepX += incX) {
    X[stepX] *= alpha;
  }
}
