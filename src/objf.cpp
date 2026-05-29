#include <complex>

#include "../include/objf.h"
#include "../include/dft2d.h"
#include "../include/math_utils.h"

typedef std::complex<float> complex;

static float* get_c(
  float* u_s, float* u_o,
  float dt, float dk,
  int nt, int nrec
)
{
  DFTOperator T = dft_operator_2d(dt, dk, nt, nrec);
  IDFTOperator TI = idft_operator_objf_2d(dt, dk, nt, nrec);

  complex* conj_im_u_s = conjugate2d(
    computeDFT(nt, nrec, u_s, T), nt, nrec
  );
  complex* im_u_o = computeDFT(nt, nrec, u_o, T);

  complex* correlation = new complex[nt * nrec];

  for (int i = 0; i < nt; ++i) {
    for (int j = 0; j < nrec; ++j) {
      int idx = i * nrec + j;

      correlation[idx] = conj_im_u_s[idx] * im_u_o[idx];
    }
  }

  return computeIDFT(correlation, TI, nt, nrec);
}

static float* get_penalty(int nt, float t0)
{
  float* P = new float[nt];

  for (int tau = -nt/2; tau < nt/2; ++tau) {
    if(std::abs(tau) <= t0) {
      P[tau] = tau;
    } else {
      P[tau] = 0.0f;
    }
  }  

  return P;
}

static float l2_sqrt(float* A, int rows, int cols)
{

//  l2[i] = np.sum((d_obs - d_calc)**2)

  float result = 0.0f;

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
       result += A[i + cols + j]; 
    }
  }

  return result;
}

float get_correlation_objf(
  float *u_s, float *u_o,
  float dt, float dk,
  int nt, int nrec, float t0
)
{
  float result = 0.0f;

  float* c = get_c(u_s, u_o, dt, dk, nt, nrec);
  float* P = get_penalty(nt, t0);

  for (int i = 0; i < nt; ++i)
  {
    float* p_c = mat_mult_scalar<float>(c, P[i], nt, nrec);
      
    result += l2_sqrt(p_c, nt, nrec);

    delete[] p_c;
  }

  return result * 0.5f;
}
