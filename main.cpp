#include <time.h>

#include "include/IO.h"
#include "include/objf.h"
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

#define PATH_CALC "data/d_calc/Seismogram_2001nt_40nrec.bin"
#define PATH_OBS "data/d_obs/seismogram_2001nt_40nrec_((250, 120))shot_33.bin"

int main()
{
  PROFILE_BEGIN();

  int nt = 2001;
  int nrec = 40;

  float dt = 1e-3;
  float dk = 5.0f;

  float t0 = 500.0f;

  float* d_calc = read2d_fortran(PATH_CALC, nt, nrec);
  float* d_obs = read2d_fortran(PATH_OBS, nt, nrec);

  float H_cor = get_correlation_objf(d_calc, d_obs, dt, dk, nt, nrec, t0);

  printf("H_cor: %f\n", H_cor);

  PROFILE_END();

  free(d_calc);
  free(d_obs);
  
  return 0;
}
