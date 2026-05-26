#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>
#include <complex>
#include <assert.h>

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

#endif // MATH_UTILS_H
