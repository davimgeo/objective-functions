#include <stdlib.h>

#include "../../include/math_utils.h"

#include "../../include/1D/ricker.h"
#include "../../include/1D/euclidian.h"
#include "../../include/1D/cross.objf.h"
#include "../../include/1D/decon_objf.h"

float* moving_rickers_l1(
  float* wavelet,
  int nt,
  int result_size
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  float phase = 0.0f;

  for (int i = 0; i < result_size; i++)
  {
    float* ricker_phase = get_ricker(
      nt,
      30.0f,
      1e-3f,
      phase
    );

    result[i] = l1_norm(
      wavelet,
      ricker_phase,
      nt
    );

    phase += 0.7f / result_size;

    free(ricker_phase);
  }

  normalize(result, result_size);

  return result;
}

float* moving_rickers_l2(
  float* wavelet,
  int nt,
  int result_size
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  float phase = 0.0f;

  for (int i = 0; i < result_size; i++)
  {
    float* ricker_phase = get_ricker(
      nt,
      30.0f,
      1e-3f,
      phase
    );

    result[i] = l2_norm(
      wavelet,
      ricker_phase,
      nt
    );

    phase += 0.7f / result_size;

    free(ricker_phase);
  }

  normalize(result, result_size);

  return result;
}

float* moving_rickers_cross(
  float* wavelet,
  int nt,
  int result_size
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  float phase = 0.0f;

  for (int i = 0; i < result_size; i++)
  {
    float* ricker_phase = get_ricker(
      nt,
      30.0f,
      1e-3f,
      phase
    );

    result[i] = get_cross_result(
      wavelet,
      ricker_phase,
      1e-3f,
      nt,
      200
    );

    phase += 0.7f / result_size;

    free(ricker_phase);
  }

  normalize(result, result_size);

  return result;
}

float* moving_rickers_decon(
  float* wavelet,
  int nt,
  int result_size
)
{
  float* result = (float*)malloc(result_size * sizeof(float));

  float phase = 0.0f;

  for (int i = 0; i < result_size; i++)
  {
    float* ricker_phase = get_ricker(
      nt,
      30.0f,
      1e-3f,
      phase
    );

    result[i] = get_decon_result(
      wavelet,
      ricker_phase,
      1e-3f,
      nt,
      200
    );

    phase += 0.7f / result_size;

    free(ricker_phase);
  }

  normalize(result, result_size);

  return result;
}

