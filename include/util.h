#ifndef DEF_UTIL_H
#define DEF_UTIL_H
#include <stdio.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

inline int min(const int a, const int b) {
  return (a < b) ? a : b;
}

#ifndef timersub
#define timersub(a, b, result) do { (result)->tv_sec = (a)->tv_sec - (b)->tv_sec; (result)->tv_usec = (a)->tv_usec - (b)->tv_usec; if ((result)->tv_usec < 0) { --(result)->tv_sec; (result)->tv_usec += 1000000; } } while (0)
#endif//timersub

double* tmp_alloc_matrix(int m, int n, double val);

void affiche(int m, int n, double *a, int lda, FILE* flux);

void init_matrix(double *mat, int m, int n, int lda, double val);

double* alloc_matrix(int m, int n);

void init_vector(double *vec, int m, double val);

double* alloc_vector(int m);

void rnd_matrix_buff(double *v, int bottom, int up, int size, int seed);

int eq_double(double a, double b, double eps);

#endif//DEF_UTIL_H
