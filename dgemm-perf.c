#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "misc.h"
#include "dgemm.h"

#ifndef SIZE
#define SIZE 10
#endif//SIZE

static double* tmp_alloc_matrix(int m, int n, double val) {
    double *res = malloc(m * n * sizeof(double));
    for (int i = 0; i < m*n; i++) {
        res[i] = val;
    }
    return res;
}

static void rnd_matrix_buff(double *v, int bottom, int up, int size, int seed) {
    srand(seed);
    for (int i = 0; i < size; ++i) {
        v[i] = rand() % (up - bottom + 1) + bottom;
    }
}

static void fill_succesion(double *v, int m, int n) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            v[i + j*m] = i*n + j + 1;
        }
    }
}

void test_matrix_product() {
    for (int size = SIZE; size <= SIZE; ++size) {
        double *A, *B, *C;
        A = tmp_alloc_matrix(size, size, 0.0); rnd_matrix_buff(A, 1, 10, size * size, 1);
        B = tmp_alloc_matrix(size, size, 0.0); rnd_matrix_buff(B, 1, 10, size * size, 1);
        C = tmp_alloc_matrix(size, size, 0.0);
        fprintf(stdout, "A %d %d\n", size, size); affiche(size, size, A, size, stdout);
        fprintf(stdout, "B %d %d\n", size, size); affiche(size, size, B, size, stdout);
        fprintf(stdout, "C %d %d\n", size, size); affiche(size, size, C, size, stdout);
        my_dgemm(COLUMN_MAJOR, 't', 'n', size, size, size, 1.0, A, size, B, size, 0.0, C, size);
        fprintf(stdout, "C %d %d\n", size, size); affiche(size, size, C, size, stdout);
        free(A);free(B);free(C);
    }
}

void test_dgemm_perf(int start, int end, int step, int nsample) {
    int size;
    struct timeval startt, endt, deltat;
    printf("size,time\n");
    for (size = start; size < end; size = (size*(100 + step))/100) {
        double *A, *B, *C;
        A = tmp_alloc_matrix(size, size, 0.0);
        B = tmp_alloc_matrix(size, size, 0.0);
        C = tmp_alloc_matrix(size, size, 0.0);
        for (int sample = 0; sample < nsample; ++sample) {
            rnd_matrix_buff(A, 1, 10, size * size, size);
            rnd_matrix_buff(B, 1, 10, size * size, size);
            gettimeofday(&startt, NULL);
            my_dgemm(COLUMN_MAJOR, 't', 'n', size, size, size, 1., A, size, B, size, 0.0, C, size);
            gettimeofday(&endt, NULL);
            timersub(&endt, &startt, &deltat);
            printf("%d,%ld.%06ld\n", size, deltat.tv_sec, deltat.tv_usec);
        }
        free(A); free(B); free(C);
    }
}

int main(int argc, char **argv) {
    // if (argc < 5) {fprintf(stderr, "argc < 5!\n"); return EXIT_FAILURE;}
    // int start, end, step, nsample;
    // start   = atoi(argv[1]);
    // end     = atoi(argv[2]);
    // step    = atoi(argv[3]);
    // nsample = atoi(argv[4]);
    test_matrix_product();
    //test_dgemm_perf(start, end, step, nsample);
    return EXIT_SUCCESS;
}