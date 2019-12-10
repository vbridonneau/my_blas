#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include "util.h"
#include "dgemm.h"
#include "algonum.h"
#include "unistd.h"
#include <stdbool.h>

double epsilon = 0.0001;
int max_width = 500;
int max_height = 500;
int max_element = 1000;
int seed = 42;
dgemm_fct_t dgemm_ref = my_dgemm_scalaire;

double * tmp_alloc_random_matrix(int m, int n) {
    double * mat = malloc(m * n * sizeof(double));
    int i;
    for(i=0;i<m*n;i++) {
        mat[i] = (rand() % (max_element+1)) * (1-(rand()%2)*2);
    }
    return mat;
}

inline void zero_matrix(int m, int n, double * mat) {
    memset(mat, m*n*sizeof(double), 0);
}

void test_dgemm(dgemm_fct_t dgemm) {
    double * a = tmp_alloc_random_matrix(max_width, max_height);
    double * b = tmp_alloc_random_matrix(max_width, max_height);
    double * c_actual = tmp_alloc_matrix(max(max_width, max_height), max(max_width, max_height), 0);
    double * c_expected = tmp_alloc_matrix(max(max_width, max_height), max(max_width, max_height), 0);

    int m = max_width;
    int n = max_height;
    int i;

    dgemm_ref(CblasColMajor, CblasTrans, CblasNoTrans, m, n, m, 1., a, n, b, n, 0.0, c_expected, n);
    dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, n, m, 1., a, n, b, n, 0.0, c_actual, n);

    bool passed = true;
    for(i=0;i<m*n;i++) {
        if(abs(c_actual[i]-c_expected[i]) > epsilon) {
            passed = false;
            break;
        }
    }

    if (passed) {
        printf("\033[0;32mPASS\033[0m");
    }
    else {
        printf("\033[0;31mFAIL\033[0m");
    }

    printf(" square\n");

    free(a);
    free(b);
    free(c_actual);
    free(c_expected);
}

int main(int argc, char * argv[]) {
	srand(seed);

    printf("---- my_dgemm_scalaire\n");
    test_dgemm(my_dgemm_scalaire);
    printf("---- my_dgemm\n");
    test_dgemm(my_dgemm);

    return EXIT_SUCCESS;
}
