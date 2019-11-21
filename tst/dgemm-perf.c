#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "util.h"
#include "dgemm.h"
#include "algonum.h"

static double* tmp_alloc_matrix(int m, int n, double val) {
    double *res = malloc(m * n * sizeof(double));
    for (int i = 0; i < m*n; i++) {
        res[i] = val;
    }
    return res;
}

static void fill_succesion(double *v, int m, int n) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            v[i + j*m] = i*n + j + 1;
        }
    }
}

#ifndef SIZE
#define SIZE 10
#endif//SIZE

const int M = 255;
const int K = 130;
const int N = 64;

void test_matrix_product() {
  /* double *A, *B, *C; */
  /* A = tmp_alloc_matrix(M, K, 0.0); rnd_matrix_buff(A, 1, 10, M * K, 1); */
  /* B = tmp_alloc_matrix(K, N, 0.0); rnd_matrix_buff(B, 1, 10, K * N, 1); */
  /* C = tmp_alloc_matrix(M, N, 0.0); */
  /* fprintf(stdout, "A %d %d\n", M, K); affiche(M, K, A, M, stdout); */
  /* fprintf(stdout, "B %d %d\n", K, N); affiche(K, N, B, K, stdout); */
  /* fprintf(stdout, "C %d %d\n", M, N); affiche(M, N, C, M, stdout); */
  /* fprintf(stdout, "a %lf\n", 1.); */
  /* fprintf(stdout, "b %lf\n", 0.); */
  /* my_dgemm(COLUMN_MAJOR, 'n', 'n', M, N, K, 1.0, A, M, B, K, 0.0, C, M); */
  /* fprintf(stdout, "C %d %d\n", M, N); affiche(M, N, C, M, stdout); */
  /* free(A);free(B);free(C); */
  char *str[] = {"OK", "NOK"};
  // int r1 = !!testall_dgemm( my_dgemm );
  // int r2 = !!testall_dgemm( my_dgemm );
}

#ifndef timersub
#define timersub(a, b, result) do { (result)->tv_sec = (a)->tv_sec - (b)->tv_sec; (result)->tv_usec = (a)->tv_usec - (b)->tv_usec; if ((result)->tv_usec < 0) { --(result)->tv_sec; (result)->tv_usec += 1000000; } } while (0)
#endif//timersub

void test_dgemm_perf(int start, int end, int step, int nsample) {
    int size;
    struct timeval startt, endt, deltat_my, deltat_mkl;
    printf("size,myperf,function\n");
    for (size = start; size < end; size = (size*(100 + step))/100) {
        double *A, *B, *C;
        A = tmp_alloc_matrix(size, size, 0.0);
        B = tmp_alloc_matrix(size, size, 0.0);
        C = tmp_alloc_matrix(size, size, 0.0);
        for (int sample = 0; sample < nsample; ++sample) {
            rnd_matrix_buff(A, 1, 10, size * size, size);
            rnd_matrix_buff(B, 1, 10, size * size, size);
	    /* Start with our version */
            gettimeofday(&startt, NULL);
            my_dgemm_scalaire(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size, 1., A, size, B, size, 0.0, C, size);
            gettimeofday(&endt, NULL);
            timersub(&endt, &startt, &deltat_my);
	    /* Followed by MKL one */
	    gettimeofday(&startt, NULL);
	    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, size, size, size, 1., A, size, B, size, 0.0, C, size);
            gettimeofday(&endt, NULL);
            timersub(&endt, &startt, &deltat_mkl);
	    double _size   = (double)size;
	    double mytime  = (double)(1000000*deltat_my.tv_sec + deltat_my.tv_usec)*1e-6;
	    double mkltime = (double)(1000000*deltat_mkl.tv_sec + deltat_mkl.tv_usec)*1e-6;
            printf("%d,%lf,my_dgemm\n%d,%lf,cblas_dgemm\n", size, flops_dgemm(_size, _size, _size)/mytime, flops_dgemm(_size, _size, _size)/mkltime);
	    fflush(stdout);
        }
        free(A); free(B); free(C);
    }
}

int main(int argc, char **argv) {
    if (argc < 5) {fprintf(stderr, "argc < 5!\n"); return EXIT_FAILURE;}
    int start, end, step, nsample;
    start   = atoi(argv[1]);
    end     = atoi(argv[2]);
    step    = atoi(argv[3]);
    nsample = atoi(argv[4]);
    // test_matrix_product();
    test_dgemm_perf(start, end, step, nsample);
    return EXIT_SUCCESS;
}
