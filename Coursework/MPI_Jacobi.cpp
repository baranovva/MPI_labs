#include <iostream>
#include <fstream>
#include <vector>
#include "mpi.h"
#include <math.h>

#define eps 1.e-7
#define N 26

int main(int argc, char *argv[]) {
	int rank, size, bottomRank, topRank;
	double time;
	std::vector<double> W, W_new, x, y;
	std::vector<int> sendYCount, sendCount, displs;
	
	auto sizeBUFFER = N*N;
	auto BUFFER = std::vector<double>(N*N);
	auto rightFunc = [](const double x, const double y) {return 0.2 * std::log(abs(x + y));};
	auto bound = [](const double x, const double y) {return 0.5 * (x + y) * (x + y);};

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if(rank == 0) {
		double temp_double;
		x.resize(N*N);
		y.resize(N*N);
		W.resize(N * N);
		W_new.resize(N * N);
		
		std::ifstream f("mesh1.msh");
		for (auto j = 0; j < N; j++) {
			for (auto i = 0; i < N; i++) {
				f >> x[N * j + i] >> y[N * j + i] >> temp_double;
			}
		}
		f.close();
	}		

	MPI_Buffer_attach(BUFFER.data(), N*N);
	
	std::vector<int> dims(1, 0);
	std::vector<int> periods(1, 0);
	
	dims[0] = size;
	periods[0] = 0;
	
	MPI_Comm commCart;
	MPI_Cart_create(MPI_COMM_WORLD, 1, dims.data(), periods.data(), 1, &commCart);
	MPI_Comm_rank(commCart, &rank);
	MPI_Cart_shift(commCart, 0, 1, &bottomRank, &topRank);	
	
	sendYCount.resize(size);
	sendCount.resize(size);
	displs.resize(size);
	if(rank == 0) {
		auto mod = N / size;
		auto remainder = N % size;
		for(auto i = 0; i < size; i++) {
			sendYCount[i] = mod;
			if(remainder != 0) {
				sendYCount[i]++;
				remainder--;
			}
			sendCount[i] = sendYCount[i] * N;
			if(i == 0)
				displs[i] = 0;
			else
				displs[i] = sendCount[i-1] + displs[i-1] ;
		}
	}
	
	MPI_Bcast(sendCount.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(sendYCount.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	auto xBuffer = std::vector<double>((sendYCount[rank] + 2)*N, 0.0);
	auto yBuffer = std::vector<double>((sendYCount[rank] + 2)*N, 0.0);

	MPI_Scatterv (x.data(), sendCount.data(), displs.data(), MPI_DOUBLE, xBuffer.data()+ N, sendCount[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv (y.data(), sendCount.data(), displs.data(), MPI_DOUBLE, yBuffer.data()+ N, sendCount[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	auto W_newBuffer = std::vector<double>((sendYCount[rank] + 2)*N, 0.0);
	for (auto j = 1; j < sendYCount[rank] + 1; j++) {
		for (auto i = 0; i < N ; i++) {
			if((bottomRank == -1 && j == 1) || (topRank == -1 && j == sendYCount[rank]) || (i == 0 || i == N-1))
				W_newBuffer[N * j + i] = bound(xBuffer[N * j + i], yBuffer[N * j + i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Status status;
	MPI_Bsend(xBuffer.data()+N, N, MPI_DOUBLE, bottomRank, 0, commCart);		
	MPI_Recv(xBuffer.data()+(sendYCount[rank]+1)*N, N, MPI_DOUBLE, topRank, 0, commCart, &status);	
	MPI_Bsend(xBuffer.data()+sendYCount[rank]*N, N, MPI_DOUBLE, topRank, 1, commCart);		
	MPI_Recv(xBuffer.data(), N, MPI_DOUBLE, bottomRank, 1, commCart, &status);	
		
	MPI_Bsend(yBuffer.data()+N, N, MPI_DOUBLE, bottomRank, 0, commCart);		
	MPI_Recv(yBuffer.data()+(sendYCount[rank]+1)*N, N, MPI_DOUBLE, topRank, 0, commCart, &status);	
	MPI_Bsend(yBuffer.data()+sendYCount[rank]*N, N, MPI_DOUBLE, topRank, 1, commCart);		
	MPI_Recv(yBuffer.data(), N, MPI_DOUBLE, bottomRank, 1, commCart, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank == 0)	time = MPI_Wtime();
	
	auto Error = 1.0;
	auto WBuffer = std::vector<double>((sendYCount[rank] + 2)*N, 0.0);
	while(Error > eps) {
		Error = 0.0;
		
		std::copy(W_newBuffer.begin(), W_newBuffer.end(), WBuffer.begin());
	
		MPI_Bsend(WBuffer.data()+N, N, MPI_DOUBLE, bottomRank, 0, commCart);		
		MPI_Recv(WBuffer.data()+(sendYCount[rank]+1)*N, N, MPI_DOUBLE, topRank, 0, commCart, &status);	
		MPI_Bsend(WBuffer.data()+sendYCount[rank]*N, N, MPI_DOUBLE, topRank, 1, commCart);		
		MPI_Recv(WBuffer.data(), N, MPI_DOUBLE, bottomRank, 1, commCart, &status);	

		auto indStart = 1;
		auto indStop= sendYCount[rank] + 1;
		
		if (bottomRank == -1) indStart = 2;
		if (topRank == -1) indStop = sendYCount[rank];
		
		auto myError = 0.0;
		for (auto j = indStart; j < indStop; j++) {
			for (auto i = 1; i < N - 1; i++) {
				auto dx_left  = (xBuffer[N * j + i] - xBuffer[N * j + i - 1]);
				auto dx_right = (xBuffer[N * j + i + 1] - xBuffer[N * j + i]);
				auto dy_down  = (yBuffer[N * j + i] - yBuffer[N * j + i - N]);
				auto dy_up    = (yBuffer[N * j + i + N] - yBuffer[N * j + i]);
				auto dxSum = dx_left + dx_right;
				auto dySum = dy_down + dy_up;
		
				W_newBuffer[N * j + i] = (W_newBuffer[N * j + i - 1] / dx_left + WBuffer[N * j + i + 1] / dx_right) / dxSum +
					(WBuffer[N * j + i - N] / dy_down + WBuffer[N * j + i + N] / dy_up) / dySum -
					0.5 * rightFunc(xBuffer[N * j + i], yBuffer[N * j + i]);
				W_newBuffer[N * j + i] /= (1.0 / dx_left / dx_right + 1.0 / dy_down / dy_up);	
				myError = std::max(myError, abs(W_newBuffer[N * j + i] - WBuffer[N * j + i]));
			}
		}	
		MPI_Allreduce(&myError, &Error, 1, MPI_DOUBLE, MPI_MAX, commCart);
	}
	
	MPI_Buffer_detach(BUFFER.data(), &sizeBUFFER);
	MPI_Gatherv(W_newBuffer.data() + N, sendCount[rank], MPI_DOUBLE, W_new.data(),
			sendCount.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(rank == 0)	{
		std::cout << "Time: " << MPI_Wtime() - time << std::endl;
		/* std::ofstream f("results.txt");
		for (auto j = 0; j < N; j++) {
			for (auto i = 0; i < N; i++) {
				f << x[N * j + i] << "\t" << y[N * j + i] << "\t" << W_new[N * j + i] << std::endl;
			}
		}
		f.close(); */
	}
	MPI_Finalize();
	return 0;
}
