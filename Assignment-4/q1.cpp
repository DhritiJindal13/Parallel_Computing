
#include <mpi.h>
#include <iostream>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int value = 0;
    int next = (rank + 1) % size;
    int prev = (rank - 1 + size) % size;

    if (rank == 0) {
        value = 100;
        value += rank;
        cout << "Process " << rank << " starts with value " << value << endl;
        MPI_Send(&value, 1, MPI_INT, next, 0, MPI_COMM_WORLD);

        MPI_Recv(&value, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "Process " << rank << " received FINAL value: " << value << endl;
    } 
    else {
        MPI_Recv(&value, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        cout << "Process " << rank << " received " << value << endl;

        value += rank;
        cout << "Process " << rank << " added rank, now " << value << endl;

        MPI_Send(&value, 1, MPI_INT, next, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}