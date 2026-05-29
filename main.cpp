#include <time.h>

#include "include/math_utils.h"
#include "include/IO.h"
#include "include/dft2d.h"
#include "include/plot.h"

struct timespec start, end;

#define PROFILE_BEGIN() \
    clock_gettime(CLOCK_MONOTONIC, &start)

#define PROFILE_END() do { \
    clock_gettime(CLOCK_MONOTONIC, &end); \
    double elapsed = (end.tv_sec - start.tv_sec) + \
                     (end.tv_nsec - start.tv_nsec) / 1e9; \
    printf("Elapsed: %.4f seconds\n", elapsed); \
} while (0)

int main()
{
  PROFILE_BEGIN();

  const char* PATH = "data/seismogram_4001nt_113nrec.bin";

  int nt = 4001;
  int nrec = 113;

  float dt = 1e-3;
  float dk = 15.0f;

  float* seismogram = read2d_fortran(PATH, nt, nrec);

  DFTOperator T = dft_operator_2d(dt, dk, nt, nrec);

  std::complex<float>* dft_seis = computeDFT(nt, nrec, seismogram, T);
  float* idft_seis = computeIDFT(dft_seis, T, nt, nrec);

  write2d(
    "data/seismogram_idft.bin",
    idft_seis,
    sizeof(float),
    nt, nrec
  );

  write2d(
    "data/seismogram_mag.bin", 
    magnitude(dft_seis, nt, nrec),
    sizeof(float),
    nt, nrec
  );

  PROFILE_END();

  //plot2d_imag(T.Tx, nrec, nrec);
  plot2d(idft_seis, nrec, nt);

  delete[] T.Tx;
  delete[] T.Tz;
  delete[] dft_seis;
  delete[] idft_seis;
  free(seismogram);
  
  return 0;
}
