#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <math.h>
#define WITHOUT_NUMPY
#include "matplotlibcpp.h"

#define PI 4 * atan(1.0)
#define NUM_FFT 2048

namespace plt = matplotlibcpp;

std::vector<std::complex<double> > fft(std::vector<std::complex<double> > x) {
  int N = x.size();
  if (N==1) return x;

  std::vector<std::complex<double> > xe(N/2,0), xo(N/2, 0), Xjo, Xjo2;
  int i;

  for(i=0;i<N;i+=2) xe[i/2] = x[i];
  for(i=1;i<N;i+=2) xo[(i-1)/2] = x[i];

  Xjo = fft(xe);
  Xjo2 = fft(xo);
  Xjo.insert (Xjo.end(), Xjo2.begin(), Xjo2.end());

  std::map<int,std::map<int,std::complex<double> > > bt;
  const std::complex<double> J(0,1);
  for(int N=2;N<=NUM_FFT;N*=2) {
    for(int k=0;k<=N/2-1;k++) {
      bt[N][k] = exp(-2*PI*k/N*J);
    }
  }

  for(i=0;i<N/2-1;i++) {
    std::complex<double> t = Xjo[i], tw = bt[N][i];
    Xjo[i] = t + tw * Xjo[i+N/2];
    Xjo[i+N/2] = t = tw * Xjo[i+N/2];
  }

  return Xjo;
}

// power spectrum ?
void powerSpectrum(std::vector<double> frame) {
  frame.resize(NUM_FFT);
  std::vector<std::complex<double> > fc(frame.begin(), frame.end());

  std::cout << frame.size() << "\n";
  std::vector<std::complex<double> > fftc = fft(fc);

  // TODO: これをすると segmentation fault が起きるので，なんとかする
  // std::cout << pow( abs(fftc[0]), 2);

  // std::vector<double> pSpectrum;
  std::vector<double> pSpectrum;
  for(int i = 0; i<NUM_FFT/2+1; i++) {
    pSpectrum[i] = pow(abs(fftc[i]),2);
  }

  // pre-emphasis and hamming filtered spectrum
  plt::plot(aSpectrum);
  plt::show();

}

void preEmphasisAndHamming(std::vector<double> frame) {
  int N = frame.size();
  static std::vector<double> f(N);

  frame[0] *= 0.08;
  for(int i = 1; i < N; i++) {
    f[i] = frame[i] - ( 0.97 * frame[i - 1]);
    frame[i] = f[i] * (0.54 - 0.46 * cos( 2.0 * PI * i / (N - 1)));
  }

  // pre-emphasis and hamming filtered spectrum
  plt::plot(frame);
  plt::show();
}

void process(int16_t* samples, size_t N) {
  std::vector<double> frame(N);

  for(int i = 0; i < N; i++) {
    frame.at(i) = samples[i];
  }

  // // original signal
  // plt::plot(frame);
  // plt::show();

  preEmphasisAndHamming(frame);
  powerSpectrum(frame);
}

int main(int argc, char ** argv) {
  if(argc < 1) return 1;
  std::ifstream wavFp;
  wavFp.open(argv[1]);
  if(!wavFp.is_open()) return 1;

  uint16_t bufferLength = 1764;
  int16_t* buffer = new int16_t[bufferLength];
  int bufferBPS = (sizeof buffer[0]);

  wavFp.seekg(6026);
  wavFp.read((char *)buffer, bufferLength*bufferBPS);

  process(buffer, bufferLength);

  wavFp.close();
}
