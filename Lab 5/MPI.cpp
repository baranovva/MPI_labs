#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "string"
#include "fstream"

#define numRows 18000
#define numCols 18000
#define NonZeroValues 10000

void SparseMatrixInCoordinateStorage(const std::vector<std::vector<double>>& matrix, std::vector<double>& values, std::vector<int>& rowInd, std::vector<int>& colInd) {
	auto numValues = 0;
	
	for (auto i = 0; i < matrix.size(); ++i) {
		for(auto j = 0; j < matrix[i].size(); ++j) {
			if(matrix[i][j] != 0) {
				++numValues;
			}
		}
	}
	
	values.resize(numValues);
	rowInd.resize(numValues);
	colInd.resize(numValues);
	
	auto index = 0;
	for (auto i = 0; i < matrix.size(); ++i) {
		for(auto j = 0; j < matrix[i].size(); ++j) {
			if(matrix[i][j] != 0) {
				values[index] = matrix[i][j];
				rowInd[index] = i;
				colInd[index] = j;
				++index;
			}
		}
	}				
}

void RandomInitializeSparseMatrix(std::vector<std::vector<double>>& matrix) {
	double val;
	int indexI, indexJ;
	
	for(auto i = 0; i < NonZeroValues; ++i) {
		val = rand();
		if(val == 0.0) --i;
		do {
			indexI = rand() % matrix.size();
			indexJ = rand() % matrix[indexI].size();
		} while (matrix[indexI][indexJ] != 0.0);
		matrix[indexI][indexJ] = val;	
	}
}

int main(int argc, char* argv[]) {
	srand((unsigned)time(NULL));
	int rank, size;
	double time;

	std::vector<double> values, valuesBuffer;
	std::vector<int> rowInd, colInd, rowIndBuffer, colIndBuffer;
	auto vector = std::vector<double>(numCols, rand());
	auto result = std::vector<double>(numRows);
	auto resultBuffer = std::vector<double>(numRows, 0.0);
	
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	auto sendCount = std::vector<int>(size);
	auto displs = std::vector<int>(size);

	if(rank == 0) {	
		auto matrix = std::vector<std::vector<double>>(numRows, std::vector<double>(numCols));	
			
		RandomInitializeSparseMatrix(matrix);
		SparseMatrixInCoordinateStorage(matrix, values, rowInd, colInd);

		int mod = values.size() / size;
		int remainder = values.size() % size;
		for(auto i = 0; i < size; i++) {
			sendCount[i] = mod;
			if(remainder != 0) {
				++sendCount[i];
				remainder--;
			}
			if (i == 0)
				displs[i] = 0;
			else
				displs[i] = sendCount[i-1] + displs[i-1]; 		
		}
		time = MPI_Wtime();								
	}
	
	MPI_Bcast(sendCount.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
		
	rowIndBuffer.resize(sendCount[rank], 0.0);
	colIndBuffer.resize(sendCount[rank], 0.0);
	valuesBuffer.resize(sendCount[rank], 0.0);
	
	MPI_Scatterv(rowInd.data(), sendCount.data(), displs.data(), MPI_INT, rowIndBuffer.data(),
		sendCount[rank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(colInd.data(), sendCount.data(), displs.data(), MPI_INT, colIndBuffer.data(),
		sendCount[rank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(values.data(), sendCount.data(), displs.data(), MPI_DOUBLE, valuesBuffer.data(),
		sendCount[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(vector.data(), numCols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	for (auto i = 0; i < sendCount[rank]; ++i) {
		resultBuffer[rowIndBuffer[i]] += valuesBuffer[i] * vector[colIndBuffer[i]];
	}
	MPI_Reduce(resultBuffer.data(), result.data(), numRows, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(rank == 0) {
		std::cout << "Time: " << MPI_Wtime() - time << std::endl;	
	}
	MPI_Finalize();
	return 0;
}
