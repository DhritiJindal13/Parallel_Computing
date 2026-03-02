#include <mpi.h>
#include <iostream>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int N = 8;
    int local_size = N / size;

    int A[8] = {1,2,3,4,5,6,7,8};
    int B[8] = {8,7,6,5,4,3,2,1};

    int local_A[8], local_B[8];
    int local_dot = 0, global_dot = 0;

    if (rank == 0) {
        MPI_Scatter(A, local_size, MPI_INT, local_A, local_size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(B, local_size, MPI_INT, local_B, local_size, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Scatter(NULL, 0, MPI_INT, local_A, local_size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(NULL, 0, MPI_INT, local_B, local_size, MPI_INT, 0, MPI_COMM_WORLD);
    }

    for (int i = 0; i < local_size; i++)
        local_dot += local_A[i] * local_B[i];

    cout << "Process " << rank << " partial dot = " << local_dot << endl;

    MPI_Reduce(&local_dot, &global_dot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
        cout << "\nFinal Dot Product = " << global_dot << " (Correct: 120)" << endl;

    MPI_Finalize();
    return 0;
}