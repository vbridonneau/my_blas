#include <stdio.h>
#include "util.h"
#include "dger.h"
#include <stdlib.h>

const int M = 3;
const int N = 7;

void test_result() {
	double *A, *X, *Y, alpha;
    A = alloc_matrix(M, N); rnd_matrix_buff(A, 1, 10, M * N, 1);
    X = alloc_vector(M);    rnd_matrix_buff(X, 1, 10, M, 1);
    Y = alloc_vector(N);    rnd_matrix_buff(Y, 1, 10, N, 1);
	alpha = (rand() % 10) + 1;
    fprintf(stdout, "A %d %d\n", M, N); affiche(M, N, A, M, stdout);
    fprintf(stdout, "X %d %d\n", M, 1); affiche(M, 1, X, M, stdout);
    fprintf(stdout, "Y %d %d\n", N, 1); affiche(N, 1, Y, N, stdout);
    fprintf(stdout, "a %lf\n", alpha);
    fprintf(stdout, "b %lf\n", 1.);
    my_dger(M, N, alpha, X, 1, Y, 1, A, M);
    fprintf(stdout, "C %d %d\n", M, N); affiche(M, N, A, M, stdout);
    free(A);free(X);free(Y);
}

int main(int argc, char **argv) {
	test_result();
	return 0;
}
