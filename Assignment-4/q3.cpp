#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL) + rank);

    int local_max = -1, local_min = 1001;
    for (int i = 0; i < 10; i++) {
        int num = rand() % 1001;
        if (num > local_max) local_max = num;
        if (num < local_min) local_min = num;
    }

    cout << "Process " << rank << " → Local Max: " << local_max 
         << " | Local Min: " << local_min << endl;

    int max_pair[2] = {local_max, rank};
    int min_pair[2] = {local_min, rank};
    int global_max[2], global_min[2];

    MPI_Reduce(max_pair, global_max, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    MPI_Reduce(min_pair, global_min, 1, MPI_2INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "\n=== FINAL RESULT ===" << endl;
        cout << "Global Maximum = " << global_max[0] 
             << " (from Process " << global_max[1] << ")" << endl;
        cout << "Global Minimum = " << global_min[0] 
             << " (from Process " << global_min[1] << ")" << endl;
    }

    MPI_Finalize();
    return 0;
}