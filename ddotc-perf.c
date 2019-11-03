#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include "ddot.h"

static double* tmp_vec_alloc(int m, double val) {
    double *res = malloc(m * sizeof(double));
    for (int i = 0; i < m; i++) {
     res[i] = val;   
    }
    return res;
}

#ifndef TIME_POINT
#define START_POINT(start) gettimeofday(&start, NULL);
#define END_POINT(end) gettimeofday(&end, NULL);
#endif//TIME_POINT

void test_ddot_perf(int start, int end, int step, int nsample) {
    int size;
    struct timeval startt, endt, deltat;
    printf("size,res,time\n");
    for (size = start; size < end; size = (size*(100 + step))/100) {
        double *X, *Y;
        for (int sample = 0; sample < nsample; ++sample) {
            X = tmp_vec_alloc(size, (double)1./(double)size);
            Y = tmp_vec_alloc(size, (double)1./(double)size);
            gettimeofday(&startt, NULL);
            double res = my_ddot(size, X, 1, Y, 1);
            gettimeofday(&endt, NULL);
            timersub(&endt, &startt, &deltat);
            printf("%d,%lf,%ld.%06ld\n", size, res, deltat.tv_sec, deltat.tv_usec);
            free(X); free(Y);
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 5) {fprintf(stderr, "argc < 5!\n"); return EXIT_FAILURE;}
    int start, end, step, nsample;
    start   = atoi(argv[1]);
    end     = atoi(argv[2]);
    step    = atoi(argv[3]);
    nsample = atoi(argv[4]);
    test_ddot_perf(start, end, step, nsample);
    return EXIT_SUCCESS;
}