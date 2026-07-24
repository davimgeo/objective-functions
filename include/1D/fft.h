#pragma once

#include <complex>

std::complex<float>* get_fft(float* arr, int N);

float* get_ifft(std::complex<float>* arr, int N);

float* fftshift(float* x, int N);
