#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

#define DEBUG() printf("%s : %s : %d\n", __FILE__, __FUNCTION__, __LINE__);

static int comp_double(double x, double y, double eps) {
  if (fabs(x - y) <= eps) return 0;
  else if (x < y) return -1;
  else return 1;
} 

void test_alloc_and_init_vec() {
  double *vec;
  int size;
  for(size = 10; size <= 1000000; size *= 10) {
    vec = alloc_vector(size);
    init_vector(vec, size, (double)size);
    /* Check all value with a certain precision */
    double eps = 1e-10;
    for (int i = 0; i < size; ++i) {
      if (comp_double(vec[i], (double)size, eps)) {goto error;}
    }
    free(vec);
  }
  printf("Test initialization et allocation de vecteur : OK\n");
  return;
  /* Error statements */
 error:
  free(vec);
  fprintf(stderr, "Erreur : le test init vecteur n'a pas fonctionne pour size = %d\n", size);
}

void test_alloc_and_init_mat() {
  double *mat;
  int m, n;
  for(m = 10; m <= 10000; m *= 10) {
    for (n = 10; n <= 10000; n *= 10) {
      mat = alloc_matrix(m, n);
      init_matrix(mat, m, n, m, (double)(m * n));
      /* Check all value with a certain precision */
      double eps = 1e-10;
      for (int i = 0; i < m * n; ++i) {
	if (comp_double(mat[i], (double)(m * n), eps)) {goto error;}
      }
      free(mat);
    }
  }
  printf("Test initialization et allocation de matrice : OK\n");
  return;
  /* Error statements */
 error:
  free(mat);
  fprintf(stderr, "Erreur : le test init matrice n'a pas fonctionne pour m = %d, n = %d\n", m, n);
}

int main(int argc, char **argv) {
  test_alloc_and_init_mat();
  test_alloc_and_init_vec();
  return 0;
}
