#include <complex.h>

#include "1D/fft.h"
#include "utils.h"

#include "plot.h"

static float get_epsilon(float complex*arr1, float complex*arr2, int size)
{
  float max = 0.0f;
  for(int i = 0; i < size; i++) 
  {
    float temp = cabsf(arr1[i] * arr2[i]);
    if (temp > max) max = temp;
  }
  return 0.01f * max;
}

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

static float* get_d(
  float* u_s, float* u_o,
  float dt, int nt
)
{
  // IFFT( (conj(A) * B) / (conj(A) * A + eps) )
  
  float complex* C_u_o = conjugate1d(get_fft(u_o, nt), nt); 
  float complex* Im_u_o = get_fft(u_o, nt);
  float complex* Im_u_s = get_fft(u_s, nt);

  float complex* d = malloc(sizeof(float complex) * nt);

  float epsilon = get_epsilon(C_u_o, Im_u_o, nt);
  for (int i = 0; i < nt; i++)
  {
    d[i] = (C_u_o[i] * Im_u_s[i]) / ((C_u_o[i] * Im_u_o[i]) + epsilon);
  }

  return get_ifft(d, nt);
}

static int initialized = 0;

float get_decon_result(
  float *u_s, float *u_o,
  float dt, int nt, float t0
)
{
  float result = 0.0f;

  float* d = get_d(u_s, u_o, dt, nt);
  float* P = get_penalty(nt, dt, t0);

  if(!initialized)
  {
    //plot1d(P, nt);
    //plot1d(d, nt);
    initialized = 1;
  }

  for (int tau = 0; tau < nt; ++tau) 
  {
    float pc = P[tau] * d[tau];

    result += pc * pc;
  }

  free(d);
  free(P);

  return (-0.5f * result);
}
