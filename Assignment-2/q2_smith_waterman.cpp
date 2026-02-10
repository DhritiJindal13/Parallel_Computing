#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

int match(char a, char b) {
    return (a == b) ? 2 : -1;
}

int main() {
    string A = "ACGTTGAC";
    string B = "ACTTG";

    int m = A.size();
    int n = B.size();

    vector<vector<int>> H(m+1, vector<int>(n+1, 0));

    int threads[] = {1, 2, 4, 8};
    double T1 = 0;

    cout << "Threads\tTime(s)\tSpeedup\tEfficiency\n";
    

    for (int t : threads) {
        omp_set_num_threads(t);

        double start = omp_get_wtime();

        for (int d = 1; d <= m + n - 1; d++) {
            #pragma omp parallel for
            for (int i = max(1, d - n); i <= min(m, d); i++) {
                int j = d - i;
                if (j >= 1 && j <= n) {
                    H[i][j] = max(0,
                        max(H[i-1][j-1] + match(A[i-1], B[j-1]),
                        max(H[i-1][j] - 1, H[i][j-1] - 1)));
                }
            }
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
