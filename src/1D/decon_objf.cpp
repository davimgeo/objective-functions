#include <complex>

#include "../../include/1D/dft.h"
#include "../../include/math_utils.h"

#include "../../include/plot.h"

static float get_epsilon(std::complex<float>* arr, int size)
{
  float max = 0.0f;
  for(int i = 0; i < size; i++) 
  {
    if (std::abs(arr[i]) > max) 
      max = std::real(arr[i]);
  }
  return 0.001f * max;
}

static float* get_penalty(int nt, float dt, float t0)
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

static float* get_d(
  float* u_s, float* u_o,
  float dt, int nt
)
{
  // IFFT( (conj(A) * B) / (conj(A) * A + eps) )
  
  std::complex<float>* C_u_o = conjugate1d(computeDFT(nt, u_o, dt), nt); 
  std::complex<float>* Im_u_o = computeDFT(nt, u_o, dt);
  std::complex<float>* Im_u_s = computeDFT(nt, u_s, dt);

  std::complex<float>* d = new std::complex<float>[nt];

  float epsilon = get_epsilon(Im_u_o, nt);
  for (int i = 0; i < nt; i++)
  {
    d[i] = (C_u_o[i] * Im_u_s[i]) / ((C_u_o[i] * Im_u_o[i]) + epsilon);
  }

  return computeIFFT(nt, d);
}

static int initialized = 0;

float get_decon_result(
  float *u_s, float *u_o,
  float dt, int nt, float t0
)
{
  float result = 0.0f;

  float* d = get_d(u_s, u_o, dt, nt);
  float* d_shift = fftshift(d, nt);
  float* P = get_penalty(nt, dt, t0);

  if(!initialized)
  {
    //plot1d(P, nt);
    plot1d(d_shift, nt);
    initialized = 1;
  }

  for (int tau = 0; tau < nt; ++tau) 
  {
    float pc = P[tau] * d_shift[tau];

    result += pc * pc;
  }

  delete[] d;
  delete[] P;

  return (0.5f * result);
}
