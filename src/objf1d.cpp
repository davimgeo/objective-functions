#include <cmath>

#include "../include/IO.h"
#include "../include/dft1d.h"
#include "../include/math_utils.h"

#include "../include/plot.h"

float* get_c(
  float* A, float* B,
  float dt, int nt
)
{
  // IFFT(conj(A) * B)

  std::complex<float>* C_A_I = conjugate1d(computeDFT(nt, A, dt), nt); 
  std::complex<float>* B_I = computeDFT(nt, B, dt);

  // cross correlation
  std::complex<float>* cross = new std::complex<float>[nt];
  for (int i = 0; i < nt; i++) cross[i] = C_A_I[i] * B_I[i];

  return computeIFFT(nt, cross);
}

float* get_penalty(int nt, float dt, float t0)
{
  float* P = new float[nt];

  for (int i = 0; i < nt; ++i) 
  {
    float tau = (i - (float)nt/2) * dt;

    if(std::abs(tau) <= t0) {
      P[i] = tau;
    } else {
      P[i] = 0.0f;
    }
  }  

  return P;
}

static int initialized = 0;

float get_correlation_objf(
  float *A, float *B,
  float dt, int nt, float t0
)
{
  // H_cor = 0.5 * sum(P(tau) * c)
  float result = 0.0f;

  float* c = get_c(A, B, dt, nt);
  float* c_shift = fftshift(c, nt);
  float* P = get_penalty(nt, dt, t0);

  if(!initialized)
  {
    //plot1d(P, nt);
    //plot1d(c_shift, nt);
    initialized = 1;
  }

  for (int tau = 0; tau < nt; ++tau) 
  {
    float pc = P[tau] * c_shift[tau];

    result += pc * pc;
  }

  delete[] c;
  delete[] P;

  return (0.5f * result);
}

float l1_norm(float* A, float* B, int size)
{
  float result = 0.0f;
  for (int i = 0; i < size; i++) 
  {
    result += fabsf(A[i] - B[i]);
  }

  return result;
}

float l2_norm(float* A, float* B, int size)
{
  float result = 0.0f;
  for (int i = 0; i < size; i++) 
  {
    result += (A[i] - B[i]) * (A[i] - B[i]);
  }

  return sqrtf(result);
}

float* get_ricker(
    int nt,
    float fmax,
    float dt,
    float phase
)
{
  float* wavelet = (float*)malloc(sizeof(float)*nt);

  float t0 = 2.0f * M_PI / fmax;

  float fc = fmax / (3.0f*sqrtf(M_PI));

  for (int i = 0; i < nt; i++) 
  {
    float t = (i * dt) - t0 - phase; 
    float arg = M_PI * (M_PI * M_PI * fc * fc * t * t);
    wavelet[i] = (1.0f - 2.0f * arg) * expf(-arg);
  }

  return wavelet;
}

void get_obj_func(
  float* wavelet,
  int nt,
  int result_size
)
{
  float* result_1 = (float*)malloc(result_size*sizeof(float));
  float* result_2 = (float*)malloc(result_size*sizeof(float));
  float* result_3 = (float*)malloc(result_size*sizeof(float));

  float phase = 0.0f;
  for (int i = 0; i < result_size; i++)
  {
    float* ricker_phase = get_ricker(nt, 30, 1e-3, phase);
 
    result_1[i] = get_correlation_objf(wavelet, ricker_phase, 1e-3, nt, 200);
    result_2[i] = l1_norm(wavelet, ricker_phase, nt);
    result_3[i] = l2_norm(wavelet, ricker_phase, nt);

    phase += 0.7f / result_size;

    free(ricker_phase);
  }

  /* Normalize result */
  normalize(result_1, result_size);
  normalize(result_2, result_size);
  normalize(result_3, result_size);

  write1d("data/1d/h_cor.bin", result_1, sizeof(float), result_size);
  write1d("data/1d/l2.bin", result_2, sizeof(float), result_size);
  write1d("data/1d/l1.bin", result_3, sizeof(float), result_size);
}

int main()
{
  int nt = 1001;
  float fmax = 30.0f;
  float dt = 1e-3;

  float* ricker = get_ricker(nt, fmax, dt, 0.3);

  get_obj_func(ricker, nt, 200);

  free(ricker);
}
