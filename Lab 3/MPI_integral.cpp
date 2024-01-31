#include <iostream>
#include <cmath>
#include <mpi.h>

#define n 300000000

int main(int argc, char *argv[]) {
    int rank, size;
    double x;
    double step = 1.0 / (double) n;
    double thread_sum = 0.0;
    double global_sum = 0.0;
    double t1, t2;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    t1 = MPI_Wtime();

    int istart = rank + 1;
    int iend = n;
    for (int i = istart; i <= iend; i += size) {
        x = (i - 0.5) * step;
        thread_sum += 1.0 / sqrt(1.0 - x * x);
    }

    MPI_Reduce(&thread_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        t2 = MPI_Wtime();
        std::cout << "Computational time = " << t2 - t1 << std::endl;
        std::cout << "Value of the integral: " << step * global_sum << std::endl;
    }

    MPI_Finalize();
    return 0;
}