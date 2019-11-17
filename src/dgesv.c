#include "dgesv.h"
#include "algonum.h"
#include <assert.h>
#include <stdbool.h>
#include "util.h"

void my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int *ldb, int *info) {
	if ( n<0 ) {
		info = -1;
	}
	else if( nrhs < 0 ){
		info = -2;
	}
	else if( lda < MAX( 1, n ) ) {
		info = -4;
	}
	else if( ldb < MAX( 1, n ) ) {
		info = -7;
	}
	if( info != 0 ) {
		// xerbla( 'DGESV ', -info ); error handling
		assert(0);
		return;
	}

	my_dgetrf( n, n, a, lda, ipiv, info );
	if( info == 0 ) {
		my_dgetrs( 'N', n, nrhs, a, lda, ipiv, b, ldb, info );
	}
}
