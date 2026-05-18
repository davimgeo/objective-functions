#include <cmath>
#include <complex.h>

#include "../include/IO.h"
#include "../include/plot.h"

int main()
{
  const char* PATH = "data/seismogram_4001nt_113nrec.bin";

  int nt = 4001;
  int nrec = 113;

  float* seismogram = read2d_fortran(PATH, nt, nrec);

  std::complex<float>* DFT;

  plot2d(seismogram, nrec, nt);

  return 0;
}

