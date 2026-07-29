#pragma once

#include <complex.h>

float complex* get_fft(float* arr, int N);

float* get_ifft(float complex* arr, int N);

float* fftshift(float* x, int N);
