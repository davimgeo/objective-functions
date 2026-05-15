#include <cmath>
#include <complex>

#define BINARY_WRITE "wb"

#define PRINT(x) printf("%f ", x);

typedef struct 
{
  std::complex<float>* DFT;
  float* freq;
} DFT;

void write1d(const char* PATH, void* arr, size_t type, int size) 
{
  FILE* bin_data = fopen(PATH, BINARY_WRITE); 
  if (bin_data == NULL) 
  {
    printf("Could not write binary file.\n");
    exit(-1);
  }

  fwrite(arr, type, size, bin_data); 

  fclose(bin_data);   
}

float* ricker(int nt, float dt, float fmax) 
{
  float* ricker = new float[nt];
  
  float t0 = 2.0f * M_PI / fmax;
  float fc = fmax / (3.0f*sqrtf(M_PI));

  for (int i = 0; i < nt; i++) 
  {
    float t = (i * dt) - t0;
    float arg = M_PI * (M_PI * M_PI * fc * fc * t * t);
    ricker[i] = (1.0f - 2.0f * arg) * expf(-arg);
  }

  return ricker;
}

DFT getDFT(int N, float* arr, float dt)
{
  DFT d;

  d.DFT = new std::complex<float>[N];
  d.freq = new float[N];

  float f_s = 1.0f / dt;

  for(int k = 0; k < N; k++)
  {
    d.DFT[k] = 0.0f;

    d.freq[k] = k * f_s / N;

    for(int n = 0; n < N; n++)
    {
      float angle = -2.0f * M_PI * k * n / N;

      d.DFT[k] += arr[n] *
              std::exp(std::complex<float>(0.0f, angle));
    }
  }

  return d;
}

float* magnitude(std::complex<float>* arr, const int SIZE)
{
  float* arr_magnitude = new float[SIZE];

  for(int i = 0; i < SIZE; i++) {
    arr_magnitude[i] = std::abs(arr[i]);
  }

  return arr_magnitude;
}

void plot1d(float* arr, int size)
{
  const char* FILE_NAME = "temp.bin";

  {
    FILE* bin_data = fopen(FILE_NAME, BINARY_WRITE);

    if (bin_data == NULL) {
      printf("Could not open the file for writing!\n");
      exit(-1);
    }

    fwrite(arr, sizeof(float), size, bin_data);

    fclose(bin_data);
  }

  FILE* gnuplot = popen("gnuplot -persistent", "w");

  if (gnuplot == NULL) {
    printf("Could not start gnuplot!\n");
    exit(-1);
  }

  fprintf(gnuplot, "plot '%s' binary format='%%float' "
                   "array=%d with lines\n",
                   FILE_NAME, size);

  fflush(gnuplot);
  pclose(gnuplot);

  remove(FILE_NAME);
}

int main()
{
  int SIZE = 1001;

  float dt = 1e-3;
  float fmax = 30.0f;

  float* wavelet = ricker(SIZE, dt, fmax);

  DFT d = getDFT(SIZE, wavelet, dt);
   
  float* dft_magnitude = magnitude(d.DFT, SIZE);

  write1d("output/ricker.bin", wavelet, sizeof(float), SIZE);
  write1d("output/mag_ricker.bin", dft_magnitude, sizeof(float), SIZE);
  write1d("output/freq.bin", d.freq, sizeof(float), SIZE);

  delete[] wavelet;
  delete[] d.DFT;
  delete[] d.freq;
  delete[] dft_magnitude;

  return 0;
}
