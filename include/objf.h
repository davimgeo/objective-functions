#ifndef OBJF_H
#define OBJF_H

float* get_c(
  float* u_s, float* u_o,
  float dt, float dk,
  int nt, int nrec
);

float* get_penalty(int nt, float t0);

float get_correlation_objf(
  float *u_s, float *u_o,
  float dt, float dk,
  int nt, int nrec,
  float t0
);


#endif // OBJ_FUNCTION_H
