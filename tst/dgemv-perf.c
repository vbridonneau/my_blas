#include "dgemv.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>

const int ROWS = 128;
const int COLS = 32;

void test_result() {
  double *X, *M, *Y, alpha, beta;
  X     = alloc_vector(ROWS); rnd_matrix_buff(X, 1, 10, ROWS, 1);
  M     = alloc_matrix(ROWS, COLS); rnd_matrix_buff(M, 1, 10, ROWS * COLS, 2);
  Y     = alloc_vector(COLS); rnd_matrix_buff(Y, 1, 10, COLS, 3);
  alpha = (rand() % 10) + 1;
  beta  = (rand() % 10) + 1;

  printf("X %d 1\n", ROWS); affiche(ROWS, 1, X, 1, stdout);
  /* The 't' only means that the operation involves A**t but we write A (not A** t) */
  printf("At %d %d\n", ROWS, COLS); affiche(ROWS, COLS, M, ROWS, stdout);
  printf("Y %d 1\n", COLS); affiche(COLS, 1, Y, 1, stdout);
  printf("a %lf\n", alpha);
  printf("b %lf\n", beta);

  my_dgemv(COLUMN_MAJOR, 't', ROWS, COLS, alpha, M, ROWS, X, 1, beta, Y, 1);

  printf("C %d 1\n", COLS); affiche(COLS, 1, Y, 1, stdout);
}

int main(int argc, char **argv) {
  test_result();
  return 0;
}
