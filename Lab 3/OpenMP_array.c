#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define N 100000000

int main() {
    int q, i, istart, iend, nthreads, tid;
    double a, t1, t2;
    double *x;

    srand(time(NULL));
    a = (double)rand() / RAND_MAX;
    x = (double*)malloc(N * sizeof(double));
    for (i = 0; i < N; i++) {
        x[i] = (double)rand() / RAND_MAX;
    }

    t1 = omp_get_wtime();
    for (q = 0; q < 2; q++) {
        #pragma omp parallel private(tid, i, istart, iend)
        {
            tid = omp_get_thread_num();
            #pragma omp single
            {
                nthreads = omp_get_num_threads();
            }
            istart = tid * N / nthreads;
            iend = istart + N / nthreads;
            for (i = istart; i < iend; i++) {
                x[i] = x[i] + a;
            }
        }
    }
    t2 = omp_get_wtime();

    printf("computational time: %f s\n", t2 - t1);
    free(x);
    return 0;
}