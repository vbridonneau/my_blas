#include <assert.h>
#include <cblas.h>
#include "dgetrf.h"
#include "algonum.h"
#include "dscal.h"
#include "dger.h"
#include "dtrsm.h"
#include "dgemm.h"

static inline int min(const int a, const int b) {
    return (a < b) ? a : b;
}

void my_dgetrf_omp(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda ) {
  int j, jb, nb;
  double one = 1.0;
  if( m == 0 || n ==0 )
    return;
  nb = min(m, n) / 2; // ilaenv( 1, 'DGETRF', ' ', m, n, -1, -1 )

  /* Fisrt step */
  my_dgetrf( Order, m, nb, a, lda);
  if (nb < n) {
    cblas_dtrsm( Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, nb, n - nb, 1., a, lda, a + nb*lda, lda);
    if (nb < m)
      my_dgemm_scal_omp(Order, CblasNoTrans, CblasNoTrans, m - nb, n - nb, nb, -1., a + nb, lda, a + nb*lda, 1., a + nb * (1 + lda), lda);
  }

  /* Second step */
  my_dgetrf( Order, m - nb, nb, a + nb * (1 + lda), lda);
  if (2*nb < n) {
    cblas_dtrsm( Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, nb, n - 2*nb, 1., a + nb * (1 + lda), lda, a + nb + 2*nb*lda, lda);
    if (2*nb < m)
      my_dgemm_scal_omp( Order, CblasNoTrans, CblasNoTrans, m - 2*nb, n - 2*nb, nb, -1., a + 2*nb + nb*lda, lda, a + nb + 2*nb*lda, lda, 1., a + 2*nb*(1 + lda), lda);
  }
  
  /* for(j = 0; j < min( m, n ); j += nb) { */
  /*   jb = min( min( m, n ) - j, nb ); */
  /*   my_dgetf2( Order, m-j, jb, a + j + j*lda, lda); */

  /*   if( j+jb < n ) { */
  /*     my_dtrsm( Order, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, jb, n-j-jb, 1., a + j*(1 + lda), lda, a + j +(j+jb)*lda, lda ); */
  /*     if( j+jb < m ) { */
  /* 	my_dgemm( Order, CblasNoTrans, CblasNoTrans, m-j-jb,n-j-jb, jb, -1., a + j+jb + j*lda, lda, a + j + (j+jb)*lda, lda, 1., a + j + jb + (j+jb)*lda, lda ); */
  /*     } */
  /*   } */
  /* } */
}
