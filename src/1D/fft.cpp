#include <fftw3.h>
#include <complex>

std::complex<float>* get_fft(float* arr, int N)
{
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  for (int i = 0; i < N; i++) 
  {
    in[i][0] = arr[i];
    in[i][1] = 0.0;
  }

  fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(p);

  auto* result = new std::complex<float>[N];

  for (int i = 0; i < N; i++) 
  {
    result[i] = std::complex<float>(
      static_cast<float>(out[i][0]),
      static_cast<float>(out[i][1])
    );
  }

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);

  return result;
}

float* get_ifft(std::complex<float>* arr, int N)
{
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  for (int i = 0; i < N; i++) 
  {
    in[i][0] = std::real(arr[i]);
    in[i][1] = std::imag(arr[i]);
  }

  fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(p);

  float* result = (float*)malloc(N * sizeof(float));
  for (int i = 0; i < N; i++) result[i] = out[i][0] / (float)N;

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);

  return result;
}

float* fftshift(float* x, int N)
{
    float* y = new float[N];

    int shift = N / 2;

    for(int i = 0; i < N; i++)
        y[i] = x[(i + shift) % N];

    return y;
}
