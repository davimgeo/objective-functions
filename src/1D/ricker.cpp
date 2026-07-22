#include <cmath>

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
