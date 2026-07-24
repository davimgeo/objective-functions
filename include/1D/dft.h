#pragma once

#include <complex>

std::complex<float>* computeDFT(int N, float* arr, float dt);

float* computeIFFT(int N, std::complex<float>* X);

float* dftshift(float* x, int N);
