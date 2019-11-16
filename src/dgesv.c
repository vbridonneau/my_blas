#include "dgesv.h"
#include "algonum.h"
#include <assert.h>
#include <stdbool.h>

void dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int *ldb, int *info) {
	// FIXME: PAS FINIIII
	  my_dgetrf(0, m, n, a, lda );
	  my_dtrsm('l', 'u', 'N', NULL, m, n, 1.0, a, lda, b, ldb);
	  my_dtrsm('l', 'l', 'N', NULL, m, n, 1.0, a, lda, b, ldb);
}
