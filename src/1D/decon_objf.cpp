#include <complex>

#include "../../include/1D/dft.h"
#include "../../include/math_utils.h"

#define EPSILON 1e-9f

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
  float* A, float* B,
  float dt, int nt
)
{
  // A = u_s; B = u_o
  // IFFT( (conj(A) * B) / (conj(A) * A + eps) )
  
  std::complex<float>* C_u_o = conjugate1d(computeDFT(nt, B, dt), nt); 
  std::complex<float>* u_o = computeDFT(nt, B, dt);
  std::complex<float>* u_s = computeDFT(nt, A, dt);

  std::complex<float>* d = new std::complex<float>[nt];
  for (int i = 0; i < nt; i++)
  {
    d[i] = (C_u_o[i] * u_s[i]) / ((C_u_o[i] * u_o[i]) + EPSILON);
  }

  return computeIFFT(nt, d);
}

static int initialized2 = 0;

float get_decon_result(
  float *A, float *B,
  float dt, int nt, float t0
)
{
  // H_cor = 0.5 * sum(P(tau) * c)
  float result = 0.0f;

  float* d = get_d(A, B, dt, nt);
  float* d_shift = fftshift(d, nt);
  float* P = get_penalty(nt, dt, t0);

  if(!initialized2)
  {
    //plot1d(P, nt);
    //plot1d(c_shift, nt);
    initialized2 = 1;
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
