#include <iostream>
#include <vector>
#include <omp.h>
using namespace std;

int main() {
    const int N = 65536;  
    vector<double> X(N), Y(N);
    double a = 2.5;

    int threads;
    double t1, t2, seq_time;

    for (int i = 0; i < N; i++) {
        X[i] = 1.0;
        Y[i] = 2.0;
    }

 
    t1 = omp_get_wtime();
    for (int i = 0; i < N; i++) {
        X[i] = a * X[i] + Y[i];
    }
    t2 = omp_get_wtime();
    seq_time = t2 - t1;
    cout << "Sequential Time = " << seq_time << " seconds\n\n";

  
    for (threads = 2; threads <= 12; threads++) {
        double par_time, speedup;

        omp_set_num_threads(threads);

        for (int i = 0; i < N; i++) {
            X[i] = 1.0;
        }

        t1 = omp_get_wtime();
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            X[i] = a * X[i] + Y[i];
        }
        t2 = omp_get_wtime();

        par_time = t2 - t1;
        speedup = seq_time / par_time;

        cout << "Threads = " << threads 
             << "   Time = " << par_time 
             << "   Speedup = " << speedup << "\n";
    }

    return 0;
}
