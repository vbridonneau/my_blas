#include "dgetrs.h"
#include <stdbool.h>
#include "algonum.h"

void my_dgetrs(
    char trans,
    int n,
    int nrhs,
    double * a,
    int lda,
    double * ipiv,
    double * b,
    int ldb,
    int info,
) {
    info = 0;
    bool notran = (trans == 'N') || (trans == 'n');

    if( !notran && ( trans != 'T' ) && ( trans != 't' ) && ( trans != 'C' ) && ( trans != 'c' ) {
        info = -1;
    }
    else if( n < 0 ) {
        info = -2;
    }
    else if( nrhs <0 ) {
        info = -3;
    }
    else if( lda < MAX( 1, n ) ) {
        info = -5;
    }
    else if( ldb < MAX( 1, n ) ) {
        info = -8;
    }
    if( info != 0 ) {
        assert(info == 0);
        //xerbla( 'DGETRS', -info )
        return;
    }
    if( n==0 || nrhs==0 )
    return;
    if( notran ) {
        my_dlaswp( nrhs, b, ldb, 1, n, ipiv, 1 );
        my_dtrsm( 'L', 'L', 'N', 'U', n, nrhs, one, a, lda, b, ldb );
        my_dtrsm( 'L', 'U', 'N', 'N', n, nrhs, one, a, lda, b, ldb );
    }
    else {
        my_dtrsm( 'L', 'U', 'T', 'N', n, nrhs, one, a, lda, b, ldb );
        my_dtrsm( 'L', 'L', 'T', 'U', n, nrhs, one, a, lda, b, ldb );
        my_dlaswp( nrhs, b, ldb, 1, n, ipiv, -1 );
    }
    return;
}
