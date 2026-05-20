#include <fftw3.h>
#include <cmath>
#include <complex>

#include "../include/IO.h"
#include "../include/plot.h"

std::complex<float>* computeDFT(int M, int N, float* arr)
{
  /* F[k,l] = \frac{1}{MN} * 
   * \sum_{m=0}^{M-1} \sum_{n=0}^{N-1} f[m,n]e^{-j2\pi( \frac{km/M} + \frac{ln}{N} )} 
   */

  std::complex<float>* DFT = new std::complex<float>[M*N];

  for(int k = 0; k < M; k++) {
    for(int l = 0; l < N; l++) {

      DFT[k + l * M] = 0.0f;

      for(int m = 0; m < M; m++) {
        for(int n = 0; n < N; n++) {
          float angle = -2.0f * M_PI *
            (((float)k * m / M) + 
             ((float)l * n / N));

          DFT[k + l * M] += arr[m + n * M] *
                  std::exp(
                      std::complex<float>(0.0f, angle)
                  );
        }
      }

      DFT[k + l * M] /= M*N;
    }
  }

  return DFT;
}

float* magnitude(fftwf_complex* arr, int height, int width)
{
  float* mag = new float[height * width];

  for(int i = 0; i < height; i++) {
    for(int j = 0; j < width; j++) {
      int idx = i + j * height;

      float re = arr[idx][0];
      float im = arr[idx][1];

      mag[idx] = std::sqrt(re*re + im*im);
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

  fftwf_complex* in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nt * nrec);
  fftwf_complex* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nt * nrec);

  for(int i = 0; i < nt * nrec; i++) {
    in[i][0] = seismogram[i];
    in[i][1] = 0.0f;
  }

  fftwf_plan p = fftwf_plan_dft_2d(nt, nrec, in, out, FFTW_FORWARD, FFTW_MEASURE);

  fftwf_execute(p);

  for(int i = 0; i < nt * nrec; i++) {
    out[i][0] /= (nt * nrec);
    out[i][1] /= (nt * nrec);
  }

  float* mag = magnitude(out, nt, nrec);

  //plot2d(mag, nrec, nt);
 
  write2d("data/seismogram_mag.bin", mag, sizeof(float), nt, nrec);

  fftwf_destroy_plan(p);
  fftwf_free(in); fftwf_free(out);

  delete[] mag;
  delete[] seismogram;

  return 0;
}

