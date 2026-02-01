#include <iostream>
#include <omp.h>
#include <iomanip>

using namespace std;

const int N = 1000;
static double A[N][N], B[N][N], C[N][N];

int main() {
    // Initialization
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = 1.0; B[i][j] = 2.0; C[i][j] = 0.0;
        }
    }

    double t1, t2, seq_time;

    // Sequential Baseline (Needed for speedup calculation)
    t1 = omp_get_wtime();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
    t2 = omp_get_wtime();
    seq_time = t2 - t1;

    // 2D Parallel Loop
    cout << "2D Parallel:" << endl;
    for (int threads = 2; threads <= 10; threads++) {
        // Reset C
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) C[i][j] = 0.0;
        
        omp_set_num_threads(threads);
        t1 = omp_get_wtime();
        
        // collapse(2) parallelizes both i and j loops
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        t2 = omp_get_wtime();
        double par_time = t2 - t1;
        cout << "Threads=" << threads << " Time=" << par_time 
             << " Speedup=" << setprecision(2) << seq_time / par_time << endl;
    }
    return 0;
}