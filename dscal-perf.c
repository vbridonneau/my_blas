#include "dscal.h"
#include "util.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>

#ifndef SIZE
#define SIZE 128
#endif//SIZE

void test_result() {
  /* Scale column scale */
  double *X, *cpy, alpha;
  X   = alloc_vector(SIZE);
  cpy = alloc_vector(SIZE);
  rnd_matrix_buff(X, 1, 10, SIZE, 1);
  memcpy(cpy, X, SIZE * sizeof(double));
  alpha = (double)((rand() % 10) + 1);
  my_dscal(SIZE, alpha, X, 1);
  for(int i = 0; i < SIZE; ++i) {
    if (!eq_double(X[i], cpy[i] * alpha, 1e-7)) {
      printf("Error : \'dscal\' not working for column vector\n");
    }
  }
  printf("\'dscal\' working for column vector\n");
  free(X); free(cpy);

  /* Scale row scale */
  double *M;
  M   = alloc_matrix(SIZE, SIZE);
  cpy = alloc_matrix(SIZE, SIZE);
  rnd_matrix_buff(X, 1, 10, SIZE * SIZE, 1);
  memcpy(cpy, M, SIZE * SIZE * sizeof(double));
  alpha = (double)((rand() % 10) + 1);
  my_dscal(SIZE, alpha, M, SIZE);
  for(int i = 0; i < SIZE; ++i) {
    /* Only the first row of M has been scale */
    if (!eq_double(X[i], cpy[i] * ((i == 0) ? alpha : 1.), 1e-7)) {
      printf("Error : \'dscal\' not working for row vector\n");
    }
  }
  printf("\'dscal\' working for row vector\n");
  free(M); free(cpy);
}

int main(int argc, char **argv) {
  test_result();
  return 0;
}
