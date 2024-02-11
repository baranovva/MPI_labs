#include <cstdlib>
#include <mpi.h>
#include <vector>

#define N 1024

MPI_Status status; 
MPI_Request request;

using std::vector;

template <class T>
class Matrix {
public:
  Matrix(int numrows, int numcols)
    :Nrow(numrows), Ncol(numcols), elements(Nrow*Ncol) {}
  Matrix(int numrows, int numcols, T* data)
    :Nrow(numrows), Ncol(numcols), elements(data, data+numrows*numcols) {}

  int rows() {return Nrow;}
  int cols() {return Ncol;}
  T operator() (int row, int col) const {return elements[Ncol*row + col];}
  T& operator() (int row, int col) {return elements[Ncol*row + col];}
  T* data() {return elements.data();}
  const vector<T>& elem() {return elements;}

private:
  int Nrow, Ncol;
  vector<T> elements;
}; 

void RandomMatrix(Matrix<double> &A, Matrix<double> &B){
  srand( time( NULL ));    
  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      A(i,j) = rand();      
    }
  }
  for (int i = 0; i < B.rows(); i++) {
    for (int j = 0; j < B.cols(); j++) {
      B(i,j) = rand();      
    }
  }
}

void VerMatrix(Matrix<double> &A, Matrix<double> &B){ 
  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      A(i,j) = 1.0;      
    }
  }
  for (int i = 0; i < B.rows(); i++) {
    for (int j = 0; j < B.cols(); j++) {
      B(i,j) = 1.0;      
    }
  }
}

int main(int argc, char *argv[]) {
  int rank, size, rowStart, rowEnd, granularity;
  double start_time;
  
  int RowsA = N;
  int ColsA = N;
  int RowsB = N;
  int ColsB = 1;  

  Matrix<double> A = Matrix<double>(RowsA, ColsA); 
  Matrix<double> B = Matrix<double>(RowsB, ColsB); 
  Matrix<double> C = Matrix<double>(RowsA, ColsB); 

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0) {
    // RandomMatrix(A, B);
    VerMatrix(A, B);
    start_time = MPI_Wtime();   
    for (int i = 1; i < size; i++) {    
      granularity = (RowsA / (size - 1)); 
      rowStart = (i - 1) * granularity;
      if (((i + 1) == size) && ((RowsA % (size - 1)) != 0)) {
  		rowEnd = RowsA;
      } else {
  		rowEnd = rowStart + granularity;
      }
      MPI_Isend(&rowStart, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &request);
      MPI_Isend(&rowEnd, 1, MPI_INT, i , 0, MPI_COMM_WORLD, &request);
      MPI_Isend(&A(rowStart,0), (rowEnd - rowStart) * ColsA, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &request);  
    }
  }
  MPI_Bcast(&B(0,0), RowsB*ColsB, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    MPI_Recv(&rowStart, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&rowEnd, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&A(rowStart,0), (rowEnd - rowStart) * ColsA, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);

    for (int i = rowStart; i < rowEnd; i++) {
      for (int j = 0; j < B.cols(); j++) {
  		for (int k = 0; k < B.rows(); k++) {
           C(i,j) += (A(i,k) * B(k,j));
  		}
      }
    }
    MPI_Isend(&rowStart, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
    MPI_Isend(&rowEnd, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&C(rowStart,0), (rowEnd - rowStart) * ColsB, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &request);
  }

  if (rank == 0) {
    for (int i = 1; i < size; i++) {
      MPI_Recv(&rowStart, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(&rowEnd, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&C(rowStart,0), (rowEnd - rowStart) * ColsB, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);     
    }
    std::cout <<"Time:" << MPI_Wtime() - start_time<< std::endl;   
    /* for (int i = 0; i < B.rows(); i++) {
    	for (int j = 0; j < B.cols(); j++) {
      		std::cout << C(i,j) << '\n';      
    		}
  		} */
  }
  MPI_Finalize();
  return 0;
}
