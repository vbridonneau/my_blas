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

void my_dgetrf2() {
    // TODO
}


void my_dgetrf(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda ) {
    int j, jb, nb;
    double one = 1.0;
    if( m == 0 || n ==0 )
    return;
    nb = 10; // ilaenv( 1, 'DGETRF', ' ', m, n, -1, -1 )
    if( nb<=1 || nb>=min( m, n ) ) {
        my_dgetrf2( m, n, a, lda, ipiv, info );
    }
    else {
        for(j=0; j<min( m, n ); j+=nb) {
            jb = min( min( m, n )-j+1, nb );
            my_dgetrf2( m-j+1, jb, &a[ j+ j*lda], lda, &ipiv[j]);

            if( j+jb <= n ) {
                my_dtrsm( CblasLeft, CblasLower, CblasNoTrans, CblasUnit, jb, n-j-jb+1, one, &a[ j+ j*lda], lda, &a[ j +(j+jb)*lda ], lda );
                if( j+jb<=m ) {
                    my_dgemm( CblasNoTrans, CblasNoTrans, m-j-jb+1,n-j-jb+1, jb, -one, &a[j+jb+ j*lda] ), lda,&a[ j+ (j+jb)*lda], lda, one, &a[ j+jb+ (j+jb)*lda],lda );
                }
            }
        }
    }
}
