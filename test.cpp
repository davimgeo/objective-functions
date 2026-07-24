#include <stdlib.h>

#include "include/IO.h"

#include "include/1D/ricker.h"
#include "include/1D/moving_rickers.h"

#include "include/plot.h"

int main()
{
  int nt = 1001;
  int result_size = 200;

  float fmax = 30.0f;
  float dt = 1e-3f;

  float t0 = 200;

  float mid = (nt * dt) / 2;
  for(int i = 0; i < result_size; i++)
  {
    float lag1 = 0.6f - i * ((0.6f - mid) / result_size);
    float lag2 = 0.05f + i * ((mid - 0.05f) / result_size);

    float* ricker1 = get_ricker(nt, fmax, dt, lag1);
    float* ricker2 = get_ricker(nt, fmax, dt, lag2);

    float* ricker = (float*)malloc(nt * sizeof(float));
    for(int t = 0; t < nt; t++) ricker[t] = ricker1[t] + ricker2[t];

    float* cross = moving_rickers_cross(ricker, nt, result_size, fmax, dt, t0);
    float* l1 = moving_rickers_l1(ricker, nt, result_size);
    float* l2 = moving_rickers_l2(ricker, nt, result_size);
    float* decon = moving_rickers_decon(ricker, nt, result_size, fmax, dt, t0);

    char path[256];

    snprintf(path, sizeof(path), "data/temp/h_cor_%d.bin", i);
    write1d(path, cross, sizeof(float), result_size);

    snprintf(path, sizeof(path), "data/temp/l1_%d.bin", i);
    write1d(path, l1, sizeof(float), result_size);

    snprintf(path, sizeof(path), "data/temp/l2_%d.bin", i);
    write1d(path, l2, sizeof(float), result_size);

    snprintf(path, sizeof(path), "data/temp/decon_%d.bin", i);
    write1d(path, decon, sizeof(float), result_size);

    free(ricker1);
    free(ricker2);
    free(ricker);
    free(cross);
    free(l1);
    free(l2);
    free(decon);

    printf("Progress: %g\n", ((float)i/result_size)*100.0f);
  }
}

