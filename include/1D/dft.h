#pragma once

#include <complex>

std::complex<float>* computeDFT(int N, float* arr, float dt);

float* computeIFFT(int N, std::complex<float>* X);

float* fftshift(float* x, int N);
