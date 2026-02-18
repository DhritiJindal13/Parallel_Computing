#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "correlate.h"

int main(int argc, char* argv[]) {
  
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <ny> <nx>\n";
        std::cerr << "  ny = number of rows (input vectors)\n";
        std::cerr << "  nx = number of columns (elements per vector)\n";
        return 1;
    }

    int ny = std::atoi(argv[1]);
    int nx = std::atoi(argv[2]);

    if (ny <= 0 || nx <= 0) {
        std::cerr << "Error: ny and nx must be positive integers.\n";
        return 1;
    }

    std::cout << "Matrix size: " << ny << " x " << nx << "\n";
    std::cout << "Output result size: " << ny << " x " << ny << "\n\n";


    float* data   = new float[ny * nx];
    float* result = new float[ny * ny];

    std::srand(42);  
    for (int i = 0; i < ny * nx; i++)
        data[i] = static_cast<float>(std::rand()) / RAND_MAX;

    auto t_start = std::chrono::high_resolution_clock::now();

    correlate(ny, nx, data, result);

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    std::cout << "correlate() completed in: " << elapsed << " seconds\n";

    std::cout << "\nSample results (first 5 diagonal entries = self-correlation):\n";
    int show = (ny < 5) ? ny : 5;
    for (int i = 0; i < show; i++)
        std::cout << "  result[" << i << "][" << i << "] = "
                  << result[i + i * ny] << "\n";
    std::cout << "(Self-correlation should be ~1.0)\n";

    delete[] data;
    delete[] result;

    return 0;
}