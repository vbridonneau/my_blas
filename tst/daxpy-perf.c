#include "daxpy.h"
#include "util.h"
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

static double* tmp_vec_alloc(int m, double val) {
    double *res = alloc_vector(m);
    init_vector(res, m, val);
    return res;
}

#ifndef TIME_POINT
#define START_POINT(start) gettimeofday(&start, NULL);
#define END_POINT(end) gettimeofday(&end, NULL);
#endif//TIME_POINT

#ifndef timersub
#define timersub(a, b, result) do { (result)->tv_sec = (a)->tv_sec - (b)->tv_sec; (result)->tv_usec = (a)->tv_usec - (b)->tv_usec; if ((result)->tv_usec < 0) { --(result)->tv_sec; (result)->tv_usec += 1000000; } } while (0)
#endif//timersub

#ifndef SIZE
#define SIZE 128
#endif//SIZE

static int cmp(double a, double b, double eps) {
  return fabs(a - b) <= eps;
}

void test_result() {
  double *X, *Y, *Z, a;
  srand(1);
  X = tmp_vec_alloc(SIZE, rand() % 10);
  Y = tmp_vec_alloc(SIZE, rand() % 10);
  Z = tmp_vec_alloc(SIZE, 1.);
  memcpy(Z, Y, SIZE * sizeof(double));
  a = rand() % 10 + 1;
  my_daxpy(SIZE, a, X, 1, Y, 1);
  double precision = 1e-9;
  for (int i = 0; i < SIZE; ++i) {
    if (!cmp(Z[i] + a*X[i], Y[i], precision)) {
      printf("daxpy : not correct\n");
      return;
    }
  }
  printf("daxpy : correct\n");
}

#ifndef timersub
#define timersub(a, b, result) do { (result)->tv_sec = (a)->tv_sec - (b)->tv_sec; (result)->tv_usec = (a)->tv_usec - (b)->tv_usec; if ((result)->tv_usec < 0) { --(result)->tv_sec; (result)->tv_usec += 1000000; } } while (0)
#endif//timersub

void test_daxpy_perf(int start, int end, int step, int nsample) {
    int size;
    struct timeval startt, endt, deltat;
    printf("size,time\n");
    for (size = start; size < end; size = (size*(100 + step))/100) {
        double *X, *Y, a;
        X = tmp_vec_alloc(size, 0.0);
        Y = tmp_vec_alloc(size, 0.0);
        a = rand() % 10 + 1;
        for (int sample = 0; sample < nsample; ++sample) {
            rnd_matrix_buff(X, 1, 10, size, size);
            rnd_matrix_buff(Y, 1, 10, size, size + 1);
            gettimeofday(&startt, NULL);
            my_daxpy(size, a, X, 1, Y, 1);
            gettimeofday(&endt, NULL);
            timersub(&endt, &startt, &deltat);
	    double _size = (double)size;
	    double time  = (double)(1000000*deltat.tv_sec + deltat.tv_usec)*1e-6;
            printf("%d,%lf\n", size, 2*_size/time);
        }
        free(X); free(Y);
    }
}

int main(int argc, char **argv) {
  if (argc < 5) {fprintf(stderr, "argc < 5!\n"); return EXIT_FAILURE;}
  int start, end, step, nsample;
  start   = atoi(argv[1]);
  end     = atoi(argv[2]);
  step    = atoi(argv[3]);
  nsample = atoi(argv[4]);
  test_daxpy_perf(start, end, step, nsample);
  //test_result();
  return 0;
}
