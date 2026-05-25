#include <fftw3.h>
#include <cmath>
#include <complex>
#include <assert.h>

#include "../include/IO.h"
#include "../include/plot.h"

struct DFTOperator
{
  std::complex<float>* Tx;
  std::complex<float>* Tz;
};

template<typename R, typename TA, typename TB>
R* mat_mult(
  const TA* A, const TB* B,
  int rowsA, int colsA, int rowsB, int colsB
)
{
  assert(colsA == rowsB);

  auto* result = new R[rowsA * colsB];

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

template<typename R, typename TA>
R* transpose(const TA* A, int rowsA, int colsA)
{
  auto* aT = new R[colsA * rowsA];

  for (int i = 0; i < rowsA; i++) {
    for (int j = 0; j < colsA; j++) {
      aT[j * rowsA + i] = A[i * colsA + j];
    }
  }  

  return aT;
}

DFTOperator dft_operator_2d(int M, int N)
{
  /* T = \sum_{m=0}^{M-1} 
  \sum_{n=0}^{N-1} e^{-j2\pi( \frac{km/M} + \frac{ln}{N} )} 
  */

  struct DFTOperator d;

  d.Tx = new std::complex<float>[N*N];
  d.Tz = new std::complex<float>[M*M];

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < M; j++) {
      float angle = -2.0f * M_PI * ((float)(i * j) / (float)M);

      d.Tz[i + j * M] =
        std::exp(std::complex<float>(0.0f, angle));
    }
  }

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
  // DFT = transpose(Tx * transpose(Tz * f))

  /*
  std::complex<float>* Tzf =
    mat_mult<std::complex<float>>(d.Tz, f, M, M, M, N);

  std::complex<float>* TxT =
    transpose<std::complex<float>>(d.Tx, N, N);

  std::complex<float>* DFT =
    mat_mult<std::complex<float>>(Tzf, TxT, M, N, N, N);

  delete[] Tzf;
  delete[] TxT;
  */
  
  std::complex<float>* Tzf =
      mat_mult<std::complex<float>>(d.Tz, f, M, M, M, N);

  std::complex<float>* TzfT =
      transpose<std::complex<float>>(Tzf, M, N);

  std::complex<float>* Tx_TzfT =
      mat_mult<std::complex<float>>(d.Tx, TzfT, N, N, N, M);

  std::complex<float>* DFT =
      transpose<std::complex<float>>(Tx_TzfT, N, M);

  delete[] Tzf;
  delete[] TzfT;
  delete[] Tx_TzfT;
  

  return DFT;
}

float* magnitude(std::complex<float>* arr, int height, int width)
{
  float* mag = new float[height * width];

  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      int idx = i + j * height;

      mag[idx] = arr[idx].real();
    }
  }

  return mag;
}

int main()
{
  const char* PATH = "data/seismogram_4001nt_113nrec.bin";

  int nt = 4001;
  int nrec = 113;

  float* seismogram = read2d_fortran(PATH, nt, nrec);

  DFTOperator T = dft_operator_2d(nt, nrec);

  plot2d_imag(T.Tz, nt, nt);

  std::complex<float>* dft_seis = computeDFT(nt, nrec, seismogram, T);

  float* mag = magnitude(dft_seis, nt, nrec);

  plot2d(mag, nrec, nt);
 
  write2d("data/seismogram_mag.bin", mag, sizeof(float), nt, nrec);

  delete[] mag;
  delete[] dft_seis;
  delete[] seismogram;

  return 0;
}

