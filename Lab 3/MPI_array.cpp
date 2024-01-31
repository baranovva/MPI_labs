#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>

int main(int argc, char* argv[]) {
    const int N = 100000000;
    int rank, size;
    double a, t1, t2;
    double *x, *local_x;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        x = new double[N];
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            x[i] = (double)rand() / RAND_MAX;
        }
	t1 = MPI_Wtime();
    }
    
    local_x = new double[N / size];
    MPI_Scatter(x, N / size, MPI_DOUBLE, local_x, N / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < N / size; i++) {
        local_x[i] = local_x[i] + a;
    }
    
    double* temp = nullptr;
    if (rank == 0) {
      temp = new double[N];
    }
    MPI_Gather(local_x, N / size, MPI_DOUBLE, temp, N / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      delete[] x;
      x = temp;
      t2 = MPI_Wtime();
      std::cout << "computational time: " << t2 - t1 << " s" << std::endl;
    }

    MPI_Finalize();
    return 0;
}