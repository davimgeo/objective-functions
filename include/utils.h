#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <complex>
#include <assert.h>

#define PROFILE_BEGIN()                                   \
    struct timespec start, end;                           \
    clock_gettime(CLOCK_MONOTONIC, &start)                \

#define PROFILE_END() do {                                \
    clock_gettime(CLOCK_MONOTONIC, &end);                 \
    double elapsed = (end.tv_sec - start.tv_sec) +        \
                     (end.tv_nsec - start.tv_nsec) / 1e9; \
    printf("Elapsed: %.4f seconds\n", elapsed);           \
} while (0)

template<typename R, typename TA, typename TB>
R* mat_mult(
  const TA* A, const TB* B,
  int rowsA, int colsA, int rowsB, int colsB
)
{
  assert(colsA == rowsB);

  auto* result = new R[rowsA * colsB];
  
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsB; j++) {
      R temp = R();

      for (int k = 0; k < colsA; k++) {
        temp += A[i * colsA + k] * B[k * colsB + j];
      }

      result[i * colsB + j] = temp;
    }
  }

  return result;
}

template<typename R, typename TA, typename TB>
R* mat_mult_scalar(
  const TA* A, const TB scalar,
  int rowsA, int colsA
)
{
  auto* result = new R[rowsA*colsA];
  
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsA; j++) {
      result[i * colsA + j] *= scalar;
    }
  }

  return result;
}

template<typename T>
T* transpose(const T* A, int rowsA, int colsA)
{
  auto* aT = new T[colsA * rowsA];

  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsA; j++) {
      aT[j * rowsA + i] = A[i * colsA + j];
    }
  }  

  return aT;
}

static inline float* magnitude(
  const std::complex<float>* arr, int height, int width
)
{
  float* mag = new float[height * width];

  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      int idx = i * width + j;

      mag[idx] = std::abs(arr[idx]);
    }
  }

  return mag;
}

static inline std::complex<float>* conjugate1d(
    const std::complex<float>* A, int size
)
{
  auto* result = new std::complex<float>[size];

  for (int i = 0; i < size; i++) {
    result[i] = std::conj(A[i]);
  }

  return result;
}

static inline std::complex<float>* conjugate2d(
  const std::complex<float>* A, int rowsA, int colsA
)
{
  int size = rowsA * colsA;
  auto* result = new std::complex<float>[size];

  for (int i = 0; i < size; i++) {
    result[i] = std::conj(A[i]);
  }

  return result;
}

static inline void normalize(float* result, int result_size)
{
  float min = result[0], max = result[0];
  for (int i = 1; i < result_size; i++)
  {
    if (result[i] < min) min = result[i];
    if (result[i] > max) max = result[i];
  }

  float range = max - min;
  if (range > 0.0f)
  {
    for (int i = 0; i < result_size; i++)
      result[i] = (result[i] - min) / range;
  }
}

#endif // MATH_UTILS_H
