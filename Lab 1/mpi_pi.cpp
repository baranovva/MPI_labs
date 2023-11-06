#include <iostream>
#include "mpi.h"

int main(int argc, char** argv) {
	double PI, x;
	int n, rank, numprocs;
	
	auto f = [](double x){return 4.0/(1.0 + x*x);};
	
	MPI_Init(&argc, &argv);
	
	double time = MPI_Wtime();
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	if(rank == 0) {	
		n = 300000000;
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	double H = 1.0/n;
	double Sum = 0.0;
	for (unsigned int i = rank; i < n; i += numprocs) {
		Sum += f(H * (i + 0.5));		
	}
	double Local_Pi = H * Sum;
	
	MPI_Reduce(&Local_Pi, &PI, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	time = MPI_Wtime() - time;
	if(rank == 0) {
		std::cout.precision(8);
		std::cout << "Pi = " << PI << std::endl;
		std::cout << "Time = " << time << std::endl;
	}
	
	MPI_Finalize();
	return 0;
}
