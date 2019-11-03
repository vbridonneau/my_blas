#include <stdio.h>
#include <stdlib.h>

void affiche(int m, int n, double *a, int lda, FILE* flux) {
  for (int l = 0; l < m; l++) {
    for (int c = 0; c < n; c++) {
      fprintf(flux, "%lf ", a[l + c*lda]);
    }
    fputc('\n', flux);
  }
}

void init_matrix(double *mat, int m, int n, int lda, double val) {
  /* Be careful with cache effect */
  for (int col = 0; col < n; ++col) {
    for (int line = 0; line < m; ++line) {
      mat[col*lda + line] = val;
    }
  }
}

double* alloc_matrix(int m, int n) {
  return (double*)malloc(m * n * sizeof(double));
}

void init_vector(double *vec, int m, double val) {
  for (int elt = 0; elt < m; ++elt) {
    vec[elt] = val;
  }
}

double* alloc_vector(int m) {
  return (double*)malloc(m * sizeof(double));
}
