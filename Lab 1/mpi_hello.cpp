#include <iostream>
#include "mpi.h"

int main(int argc, char ** argv){
	MPI_Init(&argc, &argv);
	int rank, size;
	
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	
	const int target = 0;
	if (rank == target){
		std::cout <<"Process " << rank << " of " << size << " says Hello, World!" << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}
