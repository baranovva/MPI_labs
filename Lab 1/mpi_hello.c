#include "mpi.h"
#include <stdio.h>

int main( int argc, char *argv[])
{
    int rank, size;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    printf("Process %d of %d says \"Hello, world\"\n", rank, size);

    MPI_Finalize();

    return 0;
}
