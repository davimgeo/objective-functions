#include <complex.h>
#include <stdlib.h>

#include "1D/fft.h"
#include "utils.h"

static int initialized = 0;

static float* get_penalty(int nt, float dt, float t0)
{
  float* P = (float*)malloc(nt * sizeof(float));

  for (int i = 0; i < nt; ++i) 
  {
    float tau = (i - (float)nt/2) * dt;

    if(fabs(tau) <= t0) {
      P[i] = tau;
    } else {
      P[i] = 0.0f;
    }
  }  

  return P;
}

static float* get_c(float* u_s, float* u_o, float dt, int nt)
{
  // IFFT(conj(A) * B)

  float complex* C_u_s = conjugate1d(get_fft(u_s, nt), nt); 
  float complex* Im_u_o = get_fft(u_o, nt);

  // cross correlation
  float complex* cross = malloc(sizeof(float complex) * nt);
  for (int i = 0; i < nt; i++) cross[i] = C_u_s[i] * Im_u_o[i];

  return get_ifft(cross, nt);
}

float get_cross_result(float *u_s, float *u_o, float dt, int nt, float t0)
{
  // H_cor = 0.5 * sum(P(tau) * c)
  float result = 0.0f;

  float* c = get_c(u_s, u_o, dt, nt);
  float* P = get_penalty(nt, dt, t0);

  if(!initialized)
  {
    //plot1d(P, nt);
    //plot1d(c_shift, nt);
    initialized = 1;
  }

  for (int tau = 0; tau < nt; ++tau) 
  {
    float pc = P[tau] * c[tau];

    result += pc * pc;
  }

  free(c); free(P);

  return (-0.5f * result);
}
