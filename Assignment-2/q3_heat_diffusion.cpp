#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

int main() {
    const int N = 500;
    const int STEPS = 500;

    vector<vector<double>> grid(N, vector<double>(N, 100.0));
    vector<vector<double>> new_grid(N, vector<double>(N, 0.0));

    int threads[] = {1, 2, 4, 8};
    double T1 = 0;

    cout << "Threads\tTime(s)\tSpeedup\tEfficiency\n";
  

    for (int t : threads) {
        omp_set_num_threads(t);

        double start = omp_get_wtime();

        for (int step = 0; step < STEPS; step++) {
            #pragma omp parallel for schedule(static)
            for (int i = 1; i < N-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    new_grid[i][j] = 0.25 * (
                        grid[i-1][j] + grid[i+1][j] +
                        grid[i][j-1] + grid[i][j+1]
                    );
                }
            }
            grid.swap(new_grid);
        }

        double end = omp_get_wtime();
        double time = end - start;

        if (t == 1) T1 = time;

        cout << t << "\t"
             << fixed << setprecision(6) << time << "\t"
             << setprecision(2) << (T1 / time) << "\t"
             << (T1 / time) / t << endl;
    }
    return 0;
}
