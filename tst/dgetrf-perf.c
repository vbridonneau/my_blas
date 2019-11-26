#include "dgetrf.h"
#include "algonum.h"
#include "util.h"
#include "dgemm.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

void test_dgetrf_perf(int start, int end, int step, int nsample) {
  int size;
  struct timeval startt, endt, deltat_my, deltat_mkl;
  double *A, *B, *C;
  A = tmp_alloc_matrix(size, size, 0.0);

  char *str[] = {"OK", "NOK"};

  for (size = start; size < end; size = (size*(100 + step))/100) {
    gettimeofday(&startt, NULL);
    my_dgetrf(CblasColMajor, size, size, A, size);
    gettimeofday(&endt, NULL);
    timersub(&endt, &startt, &deltat_my);
    /* Followed by MKL one */
    gettimeofday(&startt, NULL);
    //LAPACKE_dgetrf(CblasColMajor, size, size, A, size, NULL);
    gettimeofday(&endt, NULL);
    timersub(&endt, &startt, &deltat_mkl);
    double _size   = (double)size;
    double mytime  = (double)(1000000*deltat_my.tv_sec + deltat_my.tv_usec)*1e-6;
    double mkltime = (double)(1000000*deltat_mkl.tv_sec + deltat_mkl.tv_usec)*1e-6;
    printf("%d,%lf,my_dgemm\n%d,%lf,cblas_dgemm\n", size, flops_dgemm(_size, _size, _size)/mytime, size, flops_dgemm(_size, _size, _size)/mkltime);
    fflush(stdout);
  }
}

int main(int argc, char **argv) {
  /* if (argc < 5) {fprintf(stderr, "argc < 5!\n"); return EXIT_FAILURE;} */
  /* int start, end, step, nsample; */
  /* start   = atoi(argv[1]); */
  /* end     = atoi(argv[2]); */
  /* step    = atoi(argv[3]); */
  /* nsample = atoi(argv[4]); */
  /* test_dgetrf_perf(start, end, step, nsample); */
  printf("%d\n", testall_dgetrf( my_dgetrf_omp ));
  return EXIT_SUCCESS;
}
