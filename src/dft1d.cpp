#include <cmath>
#include <complex>

#include "../include/IO.h"
#include "../include/plot.h"

#define PRINT(x) printf("%f ", x);

typedef struct 
{
  std::complex<float>* DFT;
  float* freq;
} DFT;

float* ricker(int nt, float dt, float fmax) 
{
  float* ricker = new float[nt];
  
  float t0 = 2.0f * M_PI / fmax;
  float fc = fmax / (3.0f*sqrtf(M_PI));

  for (int i = 0; i < nt; i++) 
  {
    float t = (i * dt) - t0;
    float arg = M_PI * (M_PI * M_PI * fc * fc * t * t);
    ricker[i] = (1.0f - 2.0f * arg) * expf(-arg);
  }

  return ricker;
}

DFT computeDFT(int N, float* arr, float dt)
{
  DFT d;

  d.DFT = new std::complex<float>[N];
  d.freq = new float[N];

  float f_s = 1.0f / dt;

  for(int k = 0; k < N; k++)
  {
    d.DFT[k] = 0.0f;

    d.freq[k] = k * f_s / N;

    for(int n = 0; n < N; n++)
    {
      float angle = -2.0f * M_PI * k * n / N;

      d.DFT[k] += arr[n] *
              std::exp(std::complex<float>(0.0f, angle));
    }
  }

  return d;
}

std::complex<float>* computeIFFT(int N, std::complex<float>* X)
{
  std::complex<float>* x = new std::complex<float>[N]; 

  for(int n = 0; n < N; n++) {
    x[n] = 0.0f;

    for(int k = 0; k < N; k++) {
     
      float angle = 2.0f * M_PI * k * n / N;
      x[n] += (X[k] * std::exp(std::complex<float>(0.0f, angle))) / std::complex<float>(N, 0.0f);
    }
  }

  return x;
}

float* magnitude(std::complex<float>* arr, const int SIZE)
{
  float* arr_magnitude = new float[SIZE];

  for(int i = 0; i < SIZE; i++) {
    arr_magnitude[i] = std::abs(arr[i]);
  }

  return arr_magnitude;
}

int main()
{
  int SIZE = 1001;

  float dt = 2e-3;
  float fmax = 30.0f;

  int* toninho = (int*)malloc(sizeof(int) * 4.0);

  float* wavelet = ricker(SIZE, dt, fmax);

  DFT d = computeDFT(SIZE, wavelet, dt);
   
  float* dft_magnitude = magnitude(d.DFT, SIZE);

  std::complex<float>* ifft_wavelet = computeIFFT(SIZE, d.DFT);

  float* recovered_wavelet  = new float[SIZE];

  for(int i = 0; i < SIZE; i++) {
    recovered_wavelet[i] = ifft_wavelet[i].real();
  }

  plot1d_2(wavelet, recovered_wavelet, SIZE);

  //write1d("output/recovered_wavelet.bin", recovered_wavelet, sizeof(float), SIZE);
  //write1d("output/ricker.bin", wavelet, sizeof(float), SIZE);
  //write1d("output/mag_ricker.bin", dft_magnitude, sizeof(float), SIZE);
  //write1d("output/freq.bin", d.freq, sizeof(float), SIZE);

  delete[] ifft_wavelet;
  delete[] recovered_wavelet;
  delete[] wavelet;
  delete[] d.DFT;
  delete[] d.freq;
  delete[] dft_magnitude;

  return 0;
}
