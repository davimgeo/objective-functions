#ifndef DFT2D_H
#define DFT2D_H

#include <complex>

struct DFTOperator
{
  std::complex<float>* Tx;
  std::complex<float>* Tz;
};

DFTOperator dft_operator_2d(int M, int N);
DFTOperator dft_operator_2d(float dt, float dh, int M, int N);

std::complex<float>* computeDFT(int M, int N, float* f, const DFTOperator& d);
float* computeIDFT(std::complex<float>* F, const DFTOperator& d, int M, int N);

#endif // DFT2D_H
