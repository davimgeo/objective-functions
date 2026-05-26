#include <cmath>
#include <complex>
#include <time.h>
#include <assert.h>

#include "../include/math_utils.h"
#include "../include/IO.h"
#include "../include/plot.h"

struct DFTOperator
{
  std::complex<float>* Tx;
  std::complex<float>* Tz;
};

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

std::complex<float>* computeDFT(int M, int N, float* f, DFTOperator d)
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

float* computeIDFT(std::complex<float>* F, DFTOperator d, int M, int N)
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

struct timespec start, end;

int main()
{
  clock_gettime(CLOCK_MONOTONIC, &start);

  const char* PATH = "data/seismogram_4001nt_113nrec.bin";

  int nt = 4001;
  int nrec = 113;

  float* seismogram = read2d_fortran(PATH, nt, nrec);

  DFTOperator T = dft_operator_2d(nt, nrec);

  std::complex<float>* dft_seis = computeDFT(nt, nrec, seismogram, T);
  float* idft_seis = computeIDFT(dft_seis, T, nt, nrec);

  float* mag = magnitude(dft_seis, nt, nrec);

  plot2d_imag(dft_seis, nrec, nt);
 
  write2d("data/seismogram_idft.bin", idft_seis, sizeof(float), nt, nrec);
  write2d("data/seismogram_mag.bin", mag, sizeof(float), nt, nrec);

  delete[] mag;
  delete[] dft_seis;
  delete[] idft_seis;
  free(seismogram);

  // ********************************** //
  clock_gettime(CLOCK_MONOTONIC, &end);

  double elapsed = (end.tv_sec - start.tv_sec)
                  + (end.tv_nsec - start.tv_nsec) / 1e9;

  printf("Elapsed: %.4f seconds\n", elapsed);

  return 0;
}

