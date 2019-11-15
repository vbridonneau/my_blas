#include "dtrsm.h"
#include "definition.h"

int my_dtrsm(char side, char uplo, char transa, char *	diag, int m, int n, double alpha, double * a, int lda, double * b, int ldb) {
	assert(side == 'l');
	assert(uplo == 'u');
	int i, j, k;
	double lambda;

	if ((uplo == 'l') || (uplo == 'L')) {
		for(j=0; j<n; j++) {
			for(i=j+1; i<n; i++) {
				lambda = a[i+j*lda]/a[j+j*lda];
				for(k=j; k<n; k++) {
					a[i+k*lda] -= lambda * a[j+k*lda];
				}
				b[i] -= lambda*b[j];
			}
			b[j] /= a[j+j*lda];
		}
	}
	else if ((uplo == 'u') || (uplo == 'U')) {
		for(j=n-1; 0<j+1; j--) {
			for(i=j-1; 0<=i; i--) {
				lambda = a[i+j*lda]/a[j+j*lda];
				for(k=j; 0<=k; k--) {
					a[i+k*lda] -= lambda * a[j+k*lda];
				}
				b[i] -= lambda*b[j];
			}
			b[j] /= a[j+j*lda];
		}
	}
	else {
		assert(false);
	}
}
