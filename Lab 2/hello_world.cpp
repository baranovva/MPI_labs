#include <mpi.h>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::string message = "Hello, World!";
        // Блокирующая передача
        for (int i = 1; i < size; i++) {
            MPI_Send(message.c_str(), message.size() + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
    } else {
        // Блокирующий прием
        char received_message[100];
        MPI_Recv(received_message, 100, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "Process " << rank << " received message: " << received_message << std::endl;
    }
    
    if (rank == 0) {
        std::string message = "Hello, World!";
        // Неблокирующая передача
        for (int i = 1; i < size; i++) {
            MPI_Request request;
            MPI_Isend(message.c_str(), message.size() + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD, &request);
            // Проверка передачи данных до вызова MPI_Wait
            std::cout << "Sent message '" << message << "' to process " << i << " before MPI_Wait" << std::endl;
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }
    } else {
        // Неблокирующий прием
        char received_message_nonblock[100];
        MPI_Request request;
        MPI_Irecv(received_message_nonblock, 100, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        std::cout << "Process " << rank << " received message: " << received_message_nonblock << " after MPI_Wait" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
