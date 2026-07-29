#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

float complex* get_fft(float* arr, int N)
{
  fftwf_complex* in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
  fftwf_complex* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);

  for (int i = 0; i < N; i++) in[i] = arr[i] + 0.0f * I;

  fftwf_plan p = fftwf_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  fftwf_execute(p);

  float complex* result = malloc(sizeof(float complex) * N);
  memcpy(result, out, sizeof(fftwf_complex) * N);

  fftwf_destroy_plan(p);
  fftwf_free(in); fftwf_free(out);

  return result;
}

float* get_ifft(float complex* arr, int N)
{
  fftwf_complex* in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
  fftwf_complex* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);

  memcpy(in, arr, sizeof(fftwf_complex) * N);

  fftwf_plan p = fftwf_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftwf_execute(p);

  float* result = (float*)malloc(N * sizeof(float));
  for (int i = 0; i < N; i++) 
    result[i] = (out[i] + 0.0f * I) / (float)N;

  fftwf_destroy_plan(p);
  fftwf_free(in); fftwf_free(out);

  return result;
}

float* fftshift(float* x, int N)
{
  float* y = (float*)malloc(N * sizeof(float));

  int shift = N / 2;

  for(int i = 0; i < N; i++) 
    y[i] = x[(i + shift) % N];

  return y;
}
