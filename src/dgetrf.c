#include <assert.h>
#include "dgetrf.h"
#include "algonum.h"
#include "dscal.h"
#include "dger.h"

static inline int min(const int a, const int b) {
  return (a < b) ? a : b;
}

void my_dgetf2(const CBLAS_LAYOUT Order, int M, int N, double* a, int lda ) {
  int nstep = min(M, N);
  for (int j = 0; j < nstep; ++j) {
    if (j < M - 1) {
      my_dscal(M - j - 1, 1/a[j*(1 +  lda)], a + j*(1 + lda) + 1, 1);
    }
    if (j < nstep - 1) {
      my_dger(M - j - 1, N - j - 1, -1.,
	      a + j*(1 + lda) + 1, 1,
	      a + j + (j+1)*lda, lda,
	      a + (j + 1)*(1 + lda), lda);
    }
  }
}


void my_dgetrf(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda ) {
  int i,j,k;
  assert(m<=n);
  for(k=0;k<n;k++) {
    for(i=k+1; i<n; i++) {
      a[i+k*lda] /= a[k+k*lda];
      for(j=k+1;j<n;j++) {
	a[i+j*lda] -= a[i+k*lda]*a[k+j*lda];
      }
    }
  }
}
