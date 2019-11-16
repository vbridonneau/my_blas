#include "dgetrf.h"
#include "algonum.h"
#include <assert.h>

void my_dgetf2(const CBLAS_LAYOUT Order, int m, int n, double* a, int lda ) {
	int i,j,k;
	for(k=0;k<n;k++) {
		for(i=k+1; i<n; i++) {
			a[i+k*lda] /= a[k+k*lda];
			for(j=k+1;j<n;j++) {
				a[i+j*lda] -= a[i+k*lda]*a[k+j*lda];
			}
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
