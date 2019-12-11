#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"


double* tmp_alloc_matrix(int m, int n, double val) {
    double *res = malloc(m * n * sizeof(double));
    for (int i = 0; i < m*n; i++) {
        res[i] = val;
    }
    return res;
}


// Output a matrix mxn in a file
void affiche(int m, int n, double *a, int lda, FILE* flux) {
  for (int l = 0; l < m; l++) {
    for (int c = 0; c < n; c++) {
      fprintf(flux, "%lf ", a[l + c*lda]);
    }
    fputc('\n', flux);
  }
  fflush(flux);
}

// Init matrix
void init_matrix(double *mat, int m, int n, int lda, double val) {
  /* Be careful with cache effect */
  for (int col = 0; col < n; ++col) {
    for (int line = 0; line < m; ++line) {
      mat[col*lda + line] = val;
    }
  }
}

// Alloc memory space
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

void rnd_matrix_buff(double *v, int bottom, int up, int size, int seed) {
    srand(seed);
    for (int i = 0; i < size; ++i) {
        v[i] = rand() % (up - bottom + 1) + bottom;
    }
}

int eq_double(double a, double b, double eps) {
  return fabs(a - b) < eps;
}

int min(const int a, const int b) {
  return (a < b) ? a : b;
}

int max(const int a, const int b) {
    return (a > b) ? a : b;
}