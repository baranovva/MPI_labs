#include <iostream>
#include "mpi.h"

int main (int argc, char** argv) {
	int rank, size;
	
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	int COUNT = 1;
	unsigned int N = 2000000;
	int BUF = 0;
	
	MPI_Bcast(&BUF, COUNT, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier (MPI_COMM_WORLD);
	
	double Wtime = MPI_Wtime();
	
	for (unsigned int i = 0; i < N; i++)
		MPI_Bcast(&BUF, COUNT, MPI_INT, 0, MPI_COMM_WORLD);
	
	Wtime = MPI_Wtime() - Wtime;
	double Wtick = MPI_Wtick();
	
	std::cout << "Time_Wtime: " << Wtime << " Time_Wtick: "<< Wtick << std::endl;
	
	MPI_Finalize();
	return 0;
}
