
#include "correlate.h"
#include <cmath>
#include <omp.h>

#ifndef VERSION
#define VERSION 3
#endif
// VERSION 1: Sequential Baseline
#if VERSION == 1

void correlate(int ny, int nx, const float* data, float* result) {
    
    double* norm = new double[ny * nx];

    for (int y = 0; y < ny; y++) {

        double mean = 0.0;
        for (int x = 0; x < nx; x++)
            mean += data[x + y * nx];
        mean /= nx;

        double sq = 0.0;
        for (int x = 0; x < nx; x++) {
            norm[x + y * nx] = data[x + y * nx] - mean;
            sq += norm[x + y * nx] * norm[x + y * nx];
        }


        double inv = 1.0 / std::sqrt(sq);
        for (int x = 0; x < nx; x++)
            norm[x + y * nx] *= inv;
    }
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int x = 0; x < nx; x++)
                sum += norm[x + i * nx] * norm[x + j * nx];
            result[i + j * ny] = static_cast<float>(sum);
        }
    }

    delete[] norm;
}


// VERSION 2: OpenMP Parallel

#elif VERSION == 2

void correlate(int ny, int nx, const float* data, float* result) {
    double* norm = new double[ny * nx];

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; y++) {
        double mean = 0.0;
        for (int x = 0; x < nx; x++)
            mean += data[x + y * nx];
        mean /= nx;

        double sq = 0.0;
        for (int x = 0; x < nx; x++) {
            norm[x + y * nx] = data[x + y * nx] - mean;
            sq += norm[x + y * nx] * norm[x + y * nx];
        }

        double inv = 1.0 / std::sqrt(sq);
        for (int x = 0; x < nx; x++)
            norm[x + y * nx] *= inv;
    }

    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int x = 0; x < nx; x++)
                sum += norm[x + i * nx] * norm[x + j * nx];
            result[i + j * ny] = static_cast<float>(sum);
        }
    }

    delete[] norm;
}

// VERSION 3: Fully Optimized (default)

#else

void correlate(int ny, int nx, const float* data, float* result) {
    double* norm = new double[ny * nx];

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; y++) {
        double mean = 0.0;
        for (int x = 0; x < nx; x++)
            mean += data[x + y * nx];
        mean /= nx;

        double sq = 0.0;
        for (int x = 0; x < nx; x++) {
            norm[x + y * nx] = data[x + y * nx] - mean;
            sq += norm[x + y * nx] * norm[x + y * nx];
        }

        double inv = 1.0 / std::sqrt(sq);
        for (int x = 0; x < nx; x++)
            norm[x + y * nx] *= inv;
    }


    #pragma omp parallel for schedule(dynamic, 8)
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;

            #pragma omp simd reduction(+:sum)
            for (int x = 0; x < nx; x++)
                sum += norm[x + i * nx] * norm[x + j * nx];
            result[i + j * ny] = static_cast<float>(sum);
        }
    }

    delete[] norm;
}

#endif 