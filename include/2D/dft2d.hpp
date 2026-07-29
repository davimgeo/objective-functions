#ifndef DFT2D_H
#define DFT2D_H

#include <complex>

struct DFTOperator
{
  std::complex<float>* Tx;
  std::complex<float>* Tz;
};

struct IDFTOperator
{
  std::complex<float>* C_Tx;
  std::complex<float>* C_Tz_T;
};

DFTOperator dft_operator_2d(float dt, float dh, int M, int N);
IDFTOperator idft_operator_2d(float dt, float dh, int M, int N);
IDFTOperator idft_operator_objf_2d(float dt, float dh, int M, int N);

std::complex<float>* computeDFT(int M, int N, float* f, const DFTOperator& d);
float* computeIDFT(std::complex<float>* F, const IDFTOperator& d, int M, int N);

#endif // DFT2D_H
