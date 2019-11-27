#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "util.h"
#include "dgemm.h"
#include "algonum.h"
#include <omp.h>

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

void test_dgemm_perf(int start, int end, int step, int nsample, int number_thread) {
  int size;
  struct timeval startt, endt, deltat_my, deltat_mkl;
  for (size = start; size < end; size = (size*(100 + step))/100) {
    double *A, *B, *C;
    A = tmp_alloc_matrix(size, size, 0.0);
    B = tmp_alloc_matrix(size, size, 0.0);
    C = tmp_alloc_matrix(size, size, 0.0);
    
    for (int sample = 0; sample < nsample; ++sample) {
      double _size   = (double)size;
      rnd_matrix_buff(A, 1, 10, size * size, size);
      rnd_matrix_buff(B, 1, 10, size * size, size);

      /* blocked version */
      // gettimeofday(&startt, NULL);
      // my_dgemm_omp(CblasColMajor, CblasTrans, CblasNoTrans, size, size, size, 1., A, size, B, size, 1., C, size);
      // gettimeofday(&endt, NULL);
      // timersub(&endt, &startt, &deltat_my);
      // double mytime  = (double)(1000000*deltat_my.tv_sec + deltat_my.tv_usec)*1e-6;
      // printf("%d,%lf,%d\n", size, mytime, number_thread);
      
      /* scalar version */
      gettimeofday(&startt, NULL);
      my_dgemm_scal_openmp(CblasColMajor, CblasTrans, CblasNoTrans, size, size, size, 1., A, size, B, size, 0.0, C, size);
      gettimeofday(&endt, NULL);
      timersub(&endt, &startt, &deltat_my);
      double mytime  = (double)(1000000*deltat_my.tv_sec + deltat_my.tv_usec)*1e-6;
      printf("%d,%lf,%d\n", size, mytime, number_thread);


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

  int number_thread[] = {1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 60, 80};

  int i;
  printf("SIZE,TIME,THREADS\n");

  for (i=0; i<sizeof(number_thread); i++) {
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(number_thread[i]);
    test_dgemm_perf(start, end, step, nsample, number_thread[i]);
  }

  
  return EXIT_SUCCESS;
}
