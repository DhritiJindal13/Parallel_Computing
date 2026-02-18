#ifndef CORRELATE_H
#define CORRELATE_H

// Declares the correlate function interface as specified in the assignment.
// data[x + y*nx] holds element at row y, column x.
// result[i + j*ny] = Pearson correlation between row i and row j (j <= i).

void correlate(int ny, int nx, const float* data, float* result);

#endif // CORRELATE_H