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
bool verbose = true;
dgemm_fct_t dgemm_ref = my_dgemm_scalaire;


int total = 0;
int failed = 0;

#define VERBOSE(A) if(verbose) { A }

void fill_random_matrix(int m, int n, double * a, int lda) {
    int i, j;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            a[i+j*lda] = (rand() % (max_element+1)) * (1-(rand()%2)*2);
        }
     }
}

inline void zero_matrix(int m, int n, double * mat) {
    memset(mat, m*n*sizeof(double), 0);
}


bool l1_dist(int m, int n, double * a, int lda, double * b, int ldb) {
    int i, j;
    double norm = 0.0;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            norm += abs(a[i+j*lda] - b[i+j*ldb]);
        }
    }
    return norm;
}

void test_dgemm(dgemm_fct_t dgemm) {
    double * a = tmp_alloc_matrix(max_width, max_height, 0.0);
    double * b = tmp_alloc_matrix(max_width, max_height, 0.0);
    double * c_actual = tmp_alloc_matrix(max(max_width, max_height), max(max_width, max_height), 0);
    double * c_expected = tmp_alloc_matrix(max(max_width, max_height), max(max_width, max_height), 0);

    int i;
    int ta;
    int tb;
 
    CBLAS_TRANSPOSE trans[2] = {CblasTrans, CblasNoTrans};
    char * trans_text[2] = {"T", "N"};
  
    int m[3] = {max_width, rand()%max_height, rand()%max_height};
    int n[3] = {max_height, rand()%max_height, rand()%max_height}; 
    
    srand(seed);

    for(ta=0;ta<2;ta++) {
        for(tb=0;tb<2;tb++) {
            for(i=0;i<3;i++) {
                fill_random_matrix(n[i], m[i], a, n[i]);
                fill_random_matrix(m[i], n[i], b, m[i]);
                dgemm_ref(CblasColMajor, trans[ta], trans[tb], m[i], n[i], m[i], 1., a, n[i], b, m[i], 0.0, c_expected, n[i]);
                dgemm(CblasColMajor, trans[ta], trans[tb], m[i], n[i], m[i], 1., a, n[i], b, m[i], 0.0, c_actual, n[i]);
            
                VERBOSE(
                    bool passed = l1_dist(m[i], n[i], c_actual, m[i], c_expected, m[i]) < epsilon;
                    total++;
                    if (passed) {
                        printf("\033[1;32mPASS\033[0m");
                    }
                    else {
                        failed++;
                        printf("\033[1;31mFAIL\033[0m");
                    }
                    printf(" matrix %dx%d %s%s\n", m[i], n[i], trans_text[ta], trans_text[tb]);
                )
            }
        }
    }

    free(a);
    free(b);
    free(c_actual);
    free(c_expected);
}

int main(int argc, char * argv[]) {
    int opt;
    while((opt = getopt(argc, argv, "uqe:")) != -1) {
        switch (opt) {
            case 'q':
                verbose = false;
                break;
            case 'u':
                printf("Usage: %s -eq", argv[0]);
                return EXIT_SUCCESS;
            case 'e':
                epsilon=atof(optarg);
                printf("epsilon=%f\n", epsilon);
                break;
            default:
                fprintf(stderr, "invalid option %c\n", opt);
                return EXIT_FAILURE;
        }
    }

    srand(seed);

    printf("---- my_dgemm_scalaire\n");
    test_dgemm(my_dgemm_scalaire);
    printf("---- my_dgemm\n");
    test_dgemm(my_dgemm);

    if (failed) {
        printf("\n\033[0;31m%d of %d tests failed\033[0m\n", failed, total);
    }
    else {
        printf("\n\033[0;32mall tests passed successfully\033[0m\n");
    }

    return EXIT_SUCCESS;
}
