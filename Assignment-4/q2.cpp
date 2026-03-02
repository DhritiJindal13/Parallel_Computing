

#include <mpi.h>
#include <iostream>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int N = 100;
    int local_size = N / size;
    int local_array[25];
    int local_sum = 0, global_sum = 0;

    if (rank == 0) {
        int full[100];
        for (int i = 0; i < N; i++) full[i] = i + 1;
        MPI_Scatter(full, local_size, MPI_INT, local_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Scatter(NULL, 0, MPI_INT, local_array, local_size, MPI_INT, 0, MPI_COMM_WORLD);
    }

    for (int i = 0; i < local_size; i++)
        local_sum += local_array[i];

    cout << "Process " << rank << " local sum = " << local_sum << endl;

    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "\nGlobal Sum = " << global_sum << " (Correct: 5050)" << endl;
        cout << "Average = " << (global_sum / 100.0) << endl;
    }

    MPI_Finalize();
    return 0;
}