#include "algonum.h"

void my_dgemm_seq( CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
                   CBLAS_TRANSPOSE transB, const int M, const int N,
                   const int K, const double alpha, const double *A,
                   const int lda, const double *B, const int ldb,
                   const double beta, double *C, const int ldc )
{
  int m, n, k;

  if ( transA == CblasNoTrans ) {
    if ( transB == CblasNoTrans ) {
      for( k=0; k<K; k++ ) {
	double lbeta = k == 0 ? beta : 1.;
	for( m=0; m<M; m++ ) {
	  for( n=0; n<N; n++ ) {
	    C[ ldc * n + m ] = alpha * A[ lda * k + m ] * B[ ldb * n + k ]
	      +              beta * C[ ldc * n + m ];
	  }
	}
      }
    }
    else {
      for( k=0; k<K; k++ ) {
	double lbeta = k == 0 ? beta : 1.;
	for( m=0; m<M; m++ ) {
	  for( n=0; n<N; n++ ) {
	    C[ ldc * n + m ] = alpha * A[ lda * k + m ] * B[ ldb * k + n ]
	      +              beta * C[ ldc * n + m ];
	  }
	}
      }
    }
  }
  else {
    if ( transB == CblasNoTrans ) {
      for( k=0; k<K; k++ ) {
	double lbeta = k == 0 ? beta : 1.;
	for( m=0; m<M; m++ ) {
	  for( n=0; n<N; n++ ) {
	    C[ ldc * n + m ] = alpha * A[ lda * m + k ] * B[ ldb * n + k ]
	      +              beta * C[ ldc * n + m ];
	  }
	}
      }
    }
    else {
      for( k=0; k<K; k++ ) {
	double lbeta = k == 0 ? beta : 1.;
	for( m=0; m<M; m++ ) {
	  for( n=0; n<N; n++ ) {
	    C[ ldc * n + m ] = alpha * A[ lda * m + k ] * B[ ldb * k + n ]
	      +              beta * C[ ldc * n + m ];
	  }
	}
      }
    }
  }
}

void my_dgemm_scal_openmp( CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
                           CBLAS_TRANSPOSE transB, const int M, const int N,
                           const int K, const double alpha, const double *A,
                           const int lda, const double *B, const int ldb,
                           const double beta, double *C, const int ldc )
{
  my_dgemm_seq( layout, transA, transB, M, N, K,
		alpha, A, lda, B, ldb, beta, C, ldc );
}

void my_dgemm_bloc_openmp( CBLAS_LAYOUT layout, CBLAS_TRANSPOSE transA,
                           CBLAS_TRANSPOSE transB, const int M, const int N,
                           const int K, const double alpha, const double *A,
                           const int lda, const double *B, const int ldb,
                           const double beta, double *C, const int ldc )
{
  my_dgemm_seq( layout, transA, transB, M, N, K,
		alpha, A, lda, B, ldb, beta, C, ldc );
}

void my_dgemm_tiled_openmp( CBLAS_LAYOUT layout,
                            CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                            int M, int N, int K, int b,
                            double alpha, const double **A,
			    const double **B,
                            double beta,        double **C )
{
  int MT = (M + b - 1) / b;
  int NT = (N + b - 1) / b;
  int KT = (K + b - 1) / b;
  int m, n, k;
  int mm, nn, kk;
  double lbeta;

  if ( transA == CblasNoTrans ) {
    if ( transB == CblasNoTrans ) {
      for( k=0; k<KT; k++ ) {
	kk = k == (KT-1) ? K - k * b : b;
	lbeta = (k == 0) ? beta : 1.;

	for( m=0; m<MT; m++ ) {
	  mm = m == (MT-1) ? M - m * b : b;

	  for( n=0; n<NT; n++ ) {
	    nn = n == (NT-1) ? N - n * b : b;

	    cblas_dgemm( layout, transA, transB, mm, nn, kk,
			 alpha, A[ MT * k + m ], b,
			 B[ KT * n + k ], b,
			 lbeta, C[ MT * n + m ], b );
	  }
	}
      }
    }
    else {
      for( k=0; k<KT; k++ ) {
	kk = k == (KT-1) ? K - k * b : b;
	lbeta = (k == 0) ? beta : 1.;

	for( m=0; m<MT; m++ ) {
	  mm = m == (MT-1) ? M - m * b : b;

	  for( n=0; n<NT; n++ ) {
	    nn = n == (NT-1) ? N - n * b : b;

	    cblas_dgemm( layout, transA, transB, mm, nn, kk,
			 alpha, A[ MT * k + m ], b,
			 B[ NT * k + n ], b,
			 lbeta, C[ MT * n + m ], b );
	  }
	}
      }
    }
  }
  else {
    if ( transB == CblasNoTrans ) {
      for( k=0; k<KT; k++ ) {
	kk = k == (KT-1) ? K - k * b : b;
	lbeta = (k == 0) ? beta : 1.;

	for( m=0; m<MT; m++ ) {
	  mm = m == (MT-1) ? M - m * b : b;

	  for( n=0; n<NT; n++ ) {
	    nn = n == (NT-1) ? N - n * b : b;

	    cblas_dgemm( layout, transA, transB, mm, nn, kk,
			 alpha, A[ KT * m + k ], b,
			 B[ KT * n + k ], b,
			 lbeta, C[ MT * n + m ], b );
	  }
	}
      }
    }
    else {
      for( k=0; k<KT; k++ ) {
	kk = k == (KT-1) ? K - k * b : b;
	lbeta = (k == 0) ? beta : 1.;

	for( m=0; m<MT; m++ ) {
	  mm = m == (MT-1) ? M - m * b : b;

	  for( n=0; n<NT; n++ ) {
	    nn = n == (NT-1) ? N - n * b : b;

	    cblas_dgemm( layout, transA, transB, mm, nn, kk,
			 alpha, A[ KT * m + k ], b,
			 B[ NT * k + n ], b,
			 lbeta, C[ MT * n + m ], b );
	  }
	}
      }
    }
  }
}

/* To make sure we use the right prototype */
static dgemm_fct_t valig_mygemm __attribute__ ((unused)) = my_dgemm_seq;
