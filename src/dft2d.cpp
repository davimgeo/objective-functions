#include <cmath>
#include <complex>

#include "../include/math_utils.h"
#include "../include/dft2d.h"

DFTOperator dft_operator_2d(int M, int N)
{
  /* T = \sum_{m=0}^{M-1} 
  \sum_{n=0}^{N-1} e^{-j2\pi( \frac{km/M} + \frac{ln}{N} )} 
  */

  struct DFTOperator d;

  d.Tx = new std::complex<float>[N*N];
  d.Tz = new std::complex<float>[M*M];

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      float angle = -2.0f * M_PI * ((float)(i * j) / (float)M);

      d.Tz[i * M + j] =
        std::exp(std::complex<float>(0.0f, angle));
    }
  }

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      float angle = -2.0f * M_PI * ((float)(i * j) / (float)N);

      d.Tx[i * N + j] = 
        std::exp(std::complex<float>(0.0f, angle));
    }
  }

  return d;
}

DFTOperator dft_operator_2d(float dt, float dh, int M, int N)
{
  /* T = \sum_{m=0}^{M-1} 
  \sum_{n=0}^{N-1} e^{-j2\pi( \f  float* freqs_z = (float*)malloc(sizeof(float) * nt);
  */

  struct DFTOperator d;

  d.Tx = new std::complex<float>[N*N];
  d.Tz = new std::complex<float>[M*M];

  float df = 1.0f / (M*dt);
  float dk = 1.0f / (N*dh);

  #pragma omp parallel for schedule(static)
  for (int f = 0; f < M; f++) {
    for (int n = 0; n < M; n++) {
      float angle = -2.0f * M_PI * f*df * n*dt;

      d.Tz[f * M + n] =
        std::exp(std::complex<float>(0.0f, angle));
    }
  }

  #pragma omp parallel for schedule(static)
  for (int x = 0; x < N; x++) {
    for (int n = 0; n < N; n++) {
      float angle = -2.0f * M_PI * x*dk * n*dh;

      d.Tx[x * N + n] = 
        std::exp(std::complex<float>(0.0f, angle));
    }
  }

  return d;
}

std::complex<float>* computeDFT(int M, int N, float* f, const DFTOperator& d)
{
  // DFT = Tz * f * Tx^T
  
  std::complex<float>* Tzf =
    mat_mult<std::complex<float>>(d.Tz, f, M, M, M, N);

  std::complex<float>* TxT =
    transpose(d.Tx, N, N);

  std::complex<float>* DFT =
    mat_mult<std::complex<float>>(Tzf, TxT, M, N, N, N);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      DFT[i * N + j] /= M*N;
    }
  }

  delete[] Tzf;
  delete[] TxT;

  return DFT;
}

float* computeIDFT(std::complex<float>* F, const DFTOperator& d, int M, int N)
{
  // f = conj(Tz^T) * F * conj(Tx)

  std::complex<float>* conj_TzT = conjugate2d(
    transpose(d.Tz, M, M), M, M
  );

  std::complex<float>* conj_TzT_F = 
    mat_mult<std::complex<float>>(conj_TzT, F, M, M, M, N);

  std::complex<float>* conj_TzT_F_conj_Tx =
    mat_mult<std::complex<float>>(
      conj_TzT_F, conjugate2d(d.Tx, N, N), M, N, N, N
    );

  auto* IDFT = new float[M * N];

  for (int i = 0; i < M * N; ++i) {
    IDFT[i] = conj_TzT_F_conj_Tx[i].real();
  }

  return IDFT;
}


