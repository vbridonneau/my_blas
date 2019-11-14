#ifndef DEF_UTIL_H
#define DEF_UTIL_H
#include <stdio.h>

void affiche(int m, int n, double *a, int lda, FILE* flux);

void init_matrix(double *mat, int m, int n, int lda, double val);

double* alloc_matrix(int m, int n);

void init_vector(double *vec, int m, double val);

double* alloc_vector(int m);

void rnd_matrix_buff(double *v, int bottom, int up, int size, int seed);

int eq_double(double a, double b, double eps);

#endif//DEF_UTIL_H
