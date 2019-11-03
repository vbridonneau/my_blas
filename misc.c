#include <stdio.h>

void affiche(int m, int n, double *a, int lda, FILE* stream) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(stream, "%lf ", a[i + j*lda]);
        }
        fprintf(stream, "\n");
    }
}