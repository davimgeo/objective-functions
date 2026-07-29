#include <time.h>
#include <stdio.h>

#include "include/IO.h"
#include "include/plot.h"

#define PROFILE_BEGIN()                                   \
    struct timespec start, end;                           \
    clock_gettime(CLOCK_MONOTONIC, &start)                \

#define PROFILE_END() do {                                \
    clock_gettime(CLOCK_MONOTONIC, &end);                 \
    double elapsed = (end.tv_sec - start.tv_sec) +        \
                     (end.tv_nsec - start.tv_nsec) / 1e9; \
    printf("Elapsed: %.4f seconds\n", elapsed);           \
} while (0)

#define PATH_DCALC "data/d_calc/seismogram_2001nt_40nrec.bin"

#define BUFFER_SIZE 256

int main()
{
  PROFILE_BEGIN();

  int nt = 2001;
  int nrec = 40;

  float dt = 1e-3;
  float dk = 5.0f;

  float t0 = 15.0f;

  int nshots = 101;
  
  float* d_calc = read2d_fortran(PATH_DCALC, nt, nrec);
        
  float* H_cor = (float*)malloc(sizeof(float) * nshots);

  for (int i = 0; i < nshots; ++i)
  {
    char PATH_DOBS[BUFFER_SIZE];

    snprintf(
        PATH_DOBS,
        BUFFER_SIZE,
        "data/d_obs/seismogram_2001nt_40nrec_((205, 108))shot_%d.bin",
        i
    );

    float* d_obs = read2d_fortran(PATH_DOBS, nt, nrec);
    
    H_cor[i] = get_correlation_objf(
        d_calc, d_obs, dt, dk, nt, nrec, t0
    );
    
    printf("H_cor: %g, shot: %d\n", H_cor[i], i);

    free(d_obs);
  }
      
  PROFILE_END();

  char OUTPUT_PATH[BUFFER_SIZE];
  snprintf(
    OUTPUT_PATH,
    BUFFER_SIZE,
    "data/H_cor_size101_t0_%g.bin",
    t0
  );

  write1d(OUTPUT_PATH, H_cor, sizeof(float), nshots);

  free(d_calc);
  
  return 0;
}

