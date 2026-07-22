#pragma once

float* moving_rickers_l1(
  float* wavelet,
  int nt,
  int result_size
);

float* moving_rickers_l2(
float* wavelet,
int nt,
int result_size
);

float* moving_rickers_cross(
  float* wavelet,
  int nt,
  int result_size
);

float* moving_rickers_decon(
  float* wavelet,
  int nt,
  int result_size
);
