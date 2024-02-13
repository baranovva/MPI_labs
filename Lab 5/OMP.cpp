#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <string>
#include <fstream>

#define numRows 18000
#define numCols 18000
#define NonZeroValues 10000

void SparseMatrixInCoordinateStorage(const std::vector<std::vector<double>>& matrix, std::vector<double>& values, std::vector<int>& rowInd, std::vector<int>& colInd) {
    auto numValues = 0;

    for (auto i = 0; i < matrix.size(); ++i) {
        for (auto j = 0; j < matrix[i].size(); ++j) {
            if (matrix[i][j] != 0) {
                ++numValues;
            }
        }
    }

    values.resize(numValues);
    rowInd.resize(numValues);
    colInd.resize(numValues);

    auto index = 0;
    for (auto i = 0; i < matrix.size(); ++i) {
        for (auto j = 0; j < matrix[i].size(); ++j) {
            if (matrix[i][j] != 0) {
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

    for (auto i = 0; i < NonZeroValues; ++i) {
        val = rand();
        if (val == 0.0) --i;
        do {
            indexI = rand() % matrix.size();
            indexJ = rand() % matrix[indexI].size();
        } while (matrix[indexI][indexJ] != 0.0);
        matrix[indexI][indexJ] = val;
    }
}

int main() {
    srand((unsigned)time(NULL));
    double time;
    int i, tid, istart, iend, nthreads;

    std::vector<double> values;
    std::vector<int> rowInd, colInd, count, displs;
    auto matrix = std::vector<std::vector<double>>(numRows, std::vector<double>(numCols));
    auto vector = std::vector<double>(numCols, rand());
    auto result = std::vector<double>(numRows, 0.0);

    RandomInitializeSparseMatrix(matrix);
    SparseMatrixInCoordinateStorage(matrix, values, rowInd, colInd);

	#pragma omp parallel 
	{
	  #pragma omp single 
	  {
		nthreads = omp_get_num_threads();
		count.resize(nthreads, 0);
		displs.resize(nthreads, 0);
		int mod = values.size() / nthreads;
		int remainder = values.size() % nthreads;
		for(auto i = 0; i < nthreads; i++) {
		  count[i] = mod;
		  if(remainder != 0) {
		    ++count[i];
		    remainder--;
		  }
		  if(i == 0)
		    displs[i] = 0;
		  else
		    displs[i] = count[i-1] + displs[i-1];     
		}  
	  }
	}
	
    time = omp_get_wtime();
	#pragma omp parallel private(tid, istart, iend) 
	{
    	tid = omp_get_thread_num();
		istart = displs.at(tid);
		iend = istart + count.at(tid);
		for (auto i = istart; i < iend; ++i) {
		  result[rowInd[i]] += values[i] * vector[colInd[i]];
		}
  	}
    std::cout << "Time: " << omp_get_wtime() - time << std::endl;
    return 0;
}
