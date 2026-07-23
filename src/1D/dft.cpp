#include <cmath>
#include <complex>

std::complex<float>* computeDFT(int N, float* arr, float dt)
{
  std::complex<float>* DFT = new std::complex<float>[N];
  //d.freq = new float[N];

  //float f_s = 1.0f / dt;

  #pragma omp parallel for schedule(static)
  for(int k = 0; k < N; k++)
  {
    DFT[k] = 0.0f;

    //d.freq[k] = k * f_s / N;

    for(int n = 0; n < N; n++)
    {
      float angle = -2.0f * M_PI * k * n / N;

      DFT[k] += arr[n] *
              std::exp(std::complex<float>(0.0f, angle));
    }
  }

  return DFT;
}

float* computeIFFT(int N, std::complex<float>* X)
{
  float* x = new float[N]; 

  #pragma omp parallel for schedule(static)
  for(int n = 0; n < N; n++) 
  {
    x[n] = 0.0f;

    for(int k = 0; k < N; k++) 
    {
     
      float angle = 2.0f * M_PI * k * n / N;

      std::complex<float> arg = (X[k] * std::exp(std::complex<float>(0.0f, angle))) / 
        std::complex<float>(N, 0.0f);
      x[n] += arg.real();
    }
  }

  return x;
}

float* fftshift(float* x, int N)
{
    float* y = new float[N];

    int shift = N / 2;

    for(int i = 0; i < N; i++)
        y[i] = x[(i + shift) % N];

    return y;
}

float* magnitude(std::complex<float>* arr, const int SIZE)
{
  float* arr_magnitude = new float[SIZE];

  for(int i = 0; i < SIZE; i++) {
    arr_magnitude[i] = std::abs(arr[i]);
  }

  return arr_magnitude;
}

/*
int main()
{
  int SIZE = 1001;

  float dt = 2e-3;
  float fmax = 30.0f;

  int* toninho = (int*)malloc(sizeof(int) * 4.0);

  float* wavelet = ricker(SIZE, dt, fmax);

  std::complex<float>* d = computeDFT(SIZE, wavelet, dt);
   
  float* dft_magnitude = magnitude(d, SIZE);

  float* ifft_wavelet = computeIFFT(SIZE, d);

  //write1d("output/recovered_wavelet.bin", recovered_wavelet, sizeof(float), SIZE);
  //write1d("output/ricker.bin", wavelet, sizeof(float), SIZE);
  //write1d("output/mag_ricker.bin", dft_magnitude, sizeof(float), SIZE);
  //write1d("output/freq.bin", d.freq, sizeof(float), SIZE);

  delete[] ifft_wavelet;
  delete[] wavelet;
  delete[] d;
  delete[] dft_magnitude;

  return 0;
}
*/
