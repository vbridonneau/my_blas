#include <stdio.h>
#include <util.h>
#include <stdlib.h>
#include <string.h>
#include "dgemm.h"

void test_dgemm_perf(int start, int end, int step, int nsample, int maxb) {
  int size;
  struct timeval startt, endt, deltat_my, deltat_mkl;
  printf("size,perf,bloc-size\n");
  for (size = start; size < end; size = (size*(100 + step))/100) {
    double *A, *B, *C;
    A = tmp_alloc_matrix(size, size, 0.0);
    B = tmp_alloc_matrix(size, size, 0.0);
    C = tmp_alloc_matrix(size, size, 0.0);
    rnd_matrix_buff(A, 1, 10, size * size, size);
    rnd_matrix_buff(B, 1, 10, size * size, size);

    for (int b = 1; b <= maxb; ++b) {
      for (int sample = 0; sample < nsample; ++sample) {
	double _size   = (double)size;
	/* Tile A, B and C */
	double **tileA = lapack2tile( size, size, b, A, size );
	double **tileB = lapack2tile( size, size, b, A, size );
	double **tileC = lapack2tile( size, size, b, C, size );
	/* Start with our version */
	gettimeofday(&startt, NULL);
	my_dgemm_tile(CblasColMajor, CblasTrans, CblasNoTrans,
		      size, size, size, b,
		      1., (const double**)tileA, (const double**)tileB,
		      0.0, tileC);
	gettimeofday(&endt, NULL);
	timersub(&endt, &startt, &deltat_my);
	double mytime  = (double)(1000000*deltat_my.tv_sec + deltat_my.tv_usec)*1e-6;
	printf("%d,%lf,%d\n", size, flops_dgemm(_size, _size, _size)/mytime, b);
	fflush(stdout);
	/* Free tiles */
	tileFree(size, size, b, tileA);
	tileFree(size, size, b, tileB);
	tileFree(size, size, b, tileC);
      }
    }
    free(A); free(B); free(C);
  }
}


int main(int argc, char **argv) {
  if (argc == 2) {
    if(!strcmp("check", argv[1])) {
      /* Temporary */
      testall_dgemm_tiled( my_dgemm_tile );
    } else {
      fprintf(stderr, "%s : start end step nsample <max bloc size> | \"check\"!!\n", argv[0]);
      return EXIT_FAILURE;
    }
    exit(EXIT_SUCCESS);
  }

  if (argc < 6) {
    fprintf(stderr, "%s : start end step nsample <max bloc size> | \"check\"!!\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  int start, end, step, nsample, maxb;

  start   = atoi(argv[1]);
  end     = atoi(argv[2]);
  step    = atoi(argv[3]);
  nsample = atoi(argv[4]);
  maxb    = atoi(argv[5]);

  test_dgemm_perf( start, end, step, nsample, maxb );

  return EXIT_SUCCESS;
}
