#include <stdlib.h>

#include "include/IO.h"

#include "include/1D/ricker.h"
#include "include/1D/moving_rickers.h"

#include "include/plot.h"

int main()
{
  int nt = 1001;
  int result_size =200;

  float fmax = 30.0f;
  float dt = 1e-3f;

  float t0 = 500;

  #if 1
  float* ricker1 = get_ricker(nt, fmax, dt, 0.3f);
  float* ricker2 = get_ricker(nt, fmax, dt, 0.05f);

  float* ricker = (float*)malloc(nt * sizeof(float));
  for(int t = 0; t < nt; t++) ricker[t] = ricker1[t] + ricker2[t];

  write1d("data/ricker2.bin", ricker, sizeof(float), nt);
  #else
 
  float* ricker = get_ricker(nt, fmax, dt, 0.3f);
  write1d("data/ricker.bin", ricker, sizeof(float), nt);

  #endif

  float* cross = moving_rickers_cross(ricker, nt, result_size, fmax, dt, t0);
  float* l1 = moving_rickers_l1(ricker, nt, result_size);
  float* l2 = moving_rickers_l2(ricker, nt, result_size);
  float* decon = moving_rickers_decon(ricker, nt, result_size, fmax, dt, t0);

  write1d("data/1d/h_cor.bin", cross, sizeof(float), result_size);
  write1d("data/1d/l1.bin", l1, sizeof(float), result_size);
  write1d("data/1d/l2.bin", l2, sizeof(float), result_size);
  write1d("data/1d/decon.bin", decon, sizeof(float), result_size);

  free(ricker);
  free(cross);
  free(l1);
  free(l2);
  free(decon);
}
