#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

struct Particle {
    double x, y, z;
    double fx, fy, fz;
};

int main() {
    const int N = 800;
    vector<Particle> p(N);

    for (int i = 0; i < N; i++) {
        p[i].x = drand48();
        p[i].y = drand48();
        p[i].z = drand48();
        p[i].fx = p[i].fy = p[i].fz = 0.0;
    }

    int threads[] = {1, 2, 4, 8};
    double T1 = 0.0;

    cout << "Threads\tTime(s)\tSpeedup\tEfficiency\n";
   

    for (int t : threads) {
        omp_set_num_threads(t);
        double energy = 0.0;

        double start = omp_get_wtime();

        #pragma omp parallel for reduction(+:energy) schedule(dynamic)
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double dx = p[i].x - p[j].x;
                double dy = p[i].y - p[j].y;
                double dz = p[i].z - p[j].z;
                double r2 = dx*dx + dy*dy + dz*dz + 1e-12;

                double inv6 = 1.0 / (r2*r2*r2);
                double lj = 4 * (inv6*inv6 - inv6);
                energy += lj;
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
