#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "util.h"
#include "dgemm.h"
#include "algonum.h"

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
  /* my_dgemm(COLUMN_MAJOR, 'n', 'n', M, N, K, 1.0, A, M, B, K, 0.0, C, M); */
  /* free(A);free(B);free(C); */
  printf("scalar gemm :\n");
  testall_dgemm( my_dgemm_scalaire );
  printf("blocked gemm :\n");
  testall_dgemm( my_dgemm );
}

void test_dgemm_perf(int start, int end, int step, int nsample) {
  int size;
  struct timeval startt, endt, deltat_my, deltat_mkl;
  printf("size,perf,function\n");
  for (size = start; size < end; size = (size*(100 + step))/100) {
    double *A, *B, *C;
    A = tmp_alloc_matrix(size, size, 0.0);
    B = tmp_alloc_matrix(size, size, 0.0);
    C = tmp_alloc_matrix(size, size, 0.0);

    rnd_matrix_buff(A, 1, 10, size * size, size);
    rnd_matrix_buff(B, 1, 10, size * size, size);
    for (int sample = 0; sample < nsample; ++sample) {
      double _size   = (double)size;
      /* Start with our version */
      gettimeofday(&startt, NULL);
      my_dgemm_scalaire(CblasColMajor, CblasTrans, CblasNoTrans, size, size, size, 1., A, size, B, size, 0.0, C, size);
      gettimeofday(&endt, NULL);
      timersub(&endt, &startt, &deltat_my);
      double mytime  = (double)(1000000*deltat_my.tv_sec + deltat_my.tv_usec)*1e-6;
      printf("%d,%lf,my_dgemm_scalaire\n", size, flops_dgemm(_size, _size, _size)/mytime);
	    
      /* Then our blocked version */
      gettimeofday(&startt, NULL);
      my_dgemm_scal_openmp(CblasColMajor, CblasTrans, CblasNoTrans, size, size, size, 1., A, size, B, size, 0.0, C, size);
      gettimeofday(&endt, NULL);
      timersub(&endt, &startt, &deltat_my);
      mytime  = (double)(1000000*deltat_my.tv_sec + deltat_my.tv_usec)*1e-6;
      printf("%d,%lf,my_dgemm_scal_openmp\n", size, flops_dgemm(_size, _size, _size)/mytime);
      fflush(stdout);
    }
    free(A); free(B); free(C);
  }
}

int main(int argc, char **argv) {
  if (argc == 2) {
    if(!strcmp(argv[1], "check")) {
      test_matrix_product();
    } else {
      fprintf(stderr, "%s : start end step nsample | \"check\"!!\n", argv[0]);
      return EXIT_FAILURE;
    }
    exit(EXIT_SUCCESS);
  }
  if (argc < 5) {
    fprintf(stderr, "%s : start end step nsample | \"check\"!!\n", argv[0]);
    return EXIT_FAILURE;
  }
  int start, end, step, nsample;
  start   = atoi(argv[1]);
  end     = atoi(argv[2]);
  step    = atoi(argv[3]);
  nsample = atoi(argv[4]);
  test_dgemm_perf(start, end, step, nsample);
  test_matrix_product();
  return EXIT_SUCCESS;
}
