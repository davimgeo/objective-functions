#include <cmath>
#include <complex>

#include "../../include/1D/dft.h"
#include "../../include/math_utils.h"

static int initialized = 0;

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

static float* get_c(
  float* A, float* B,
  float dt, int nt
)
{
  // A = u_s; B = u_o
  // IFFT(conj(A) * B)

  std::complex<float>* C_u_s = conjugate1d(computeDFT(nt, A, dt), nt); 
  std::complex<float>* u_o = computeDFT(nt, B, dt);

  // cross correlation
  std::complex<float>* cross = new std::complex<float>[nt];
  for (int i = 0; i < nt; i++) cross[i] = C_u_s[i] * u_o[i];

  return computeIFFT(nt, cross);
}

float get_cross_result(
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
